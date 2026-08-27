cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

set(CMAKE_CXX_STANDARD 17)

project(cafnusyst VERSION 1.00.00 LANGUAGES CXX)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE DebWithRelInfo)
endif()

find_package(nusystematics 2.00.03 REQUIRED)

find_package(duneanaobj REQUIRED)

message(STATUS "Building against duneanaobj")

# StandardRecord and its ROOT dictionary are always required.
set(ANAOBJ_LIBS
    duneanaobj::StandardRecord
    duneanaobj::StandardRecord_dict
)

# StandardRecordProxy / StandardRecordFlat are not exported as CMake targets by
# every duneanaobj install: a default duneanaobj build (SKIP_CET=OFF, the
# default) compiles these two libraries via cetmodules' cet_make_library(),
# whose export doesn't land in the duneanaobjTargets export set that
# duneanaobj::StandardRecordProxy/StandardRecordFlat come from -- only a
# SKIP_CET=ON duneanaobj build exports them that way. Either way the compiled
# library genuinely exists on disk (Instantiations.cxx / SRProxy.cxx /
# FlatRecord.cxx are only ever compiled into these libraries), so when the
# target isn't exported, fall back to locating the library file directly by
# name, using duneanaobj_DIR (set by the find_package(duneanaobj) above) as a
# search hint. This matters because some of the missing symbols (e.g.
# caf::GetCAFType) are plain compiled functions, not templates -- headers
# alone can't supply them.
#
# Kept out of ANAOBJ_LIBS deliberately: on a default duneanaobj build these
# resolve to *static* archives (.a). A static archive only contributes the
# object files needed to resolve symbols already referenced by whatever it's
# linked into -- and it gets a private copy of that object code embedded
# there. Since that object code includes ROOT dictionary-registration statics
# (TGenericClassInfo, one per caf::Proxy<T> specialization), multiple copies
# loaded into the same process register/unregister the same TClass entries
# independently -- corrupting ROOT's bookkeeping and crashing with heap
# corruption (malloc_consolidate) when the second copy's static destructors
# run at exit. So each archive must be linked into exactly the one library
# that actually references its symbols directly (StandardRecordProxy into
# cafnusyst_utility, via NuTreeHelper's caf::Proxy<T> use; StandardRecordFlat
# into cafnusyst_interface, via WeightUpdater's caf::FlatStandardRecord use --
# see each target's own CMakeLists.txt PRIVATE section), and picked up
# transitively everywhere else via the existing cafnusyst::utility /
# cafnusyst::interface dependencies. Linking a given archive into more than
# one of those two places, or into a library that doesn't itself reference
# its symbols, either reintroduces the duplicate-registration crash above or
# silently links in nothing at all (the archive contributes no object files
# if nothing in that link references them).
function(cafnusyst_resolve_anaobj_lib _libname _outvar)
  if(TARGET duneanaobj::${_libname})
    set(${_outvar} duneanaobj::${_libname} PARENT_SCOPE)
  else()
    find_library(${_libname}_LIBRARY
      NAMES ${_libname}
      HINTS ${duneanaobj_DIR}/../../..
      PATH_SUFFIXES lib lib64)
    if(${_libname}_LIBRARY)
      message(STATUS "duneanaobj target duneanaobj::${_libname} not exported; "
                      "linking ${${_libname}_LIBRARY} directly")
      set(${_outvar} ${${_libname}_LIBRARY} PARENT_SCOPE)
    else()
      message(STATUS "duneanaobj target/library ${_libname} not found; "
                      "relying on headers/explicit instantiation")
      set(${_outvar} "" PARENT_SCOPE)
    endif()
  endif()
endfunction()

cafnusyst_resolve_anaobj_lib(duneanaobj_StandardRecordProxy ANAOBJ_PROXY_LIB)
cafnusyst_resolve_anaobj_lib(duneanaobj_StandardRecordFlat  ANAOBJ_FLAT_LIB)

find_package(GENIE3 REQUIRED)
if(NOT TARGET GENIE3::All)
  message(FATAL_ERROR "Expected find_package(GENIE3 REQUIRED) call to set up target GENIE3::All.")
endif()

###### Compiler set up
add_library(cafnusyst_dependencies INTERFACE)
target_include_directories(cafnusyst_dependencies INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include> )
set_target_properties(cafnusyst_dependencies PROPERTIES EXPORT_NAME commondeps)
add_library(cafnusyst::commondeps ALIAS cafnusyst_dependencies)

install(TARGETS cafnusyst_dependencies EXPORT cafnusyst-targets)

# Preserve the RPATH for external libraries (ROOT, GENIE, yaml-cpp) during install
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Ensure installed binaries can find internal shared libraries installed to lib/
# Note: Use CMAKE_INSTALL_LIBDIR since Eigen3 sets it, or hardcode "lib"
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# Also RPATH the compiler's own implicit runtime libs (e.g. libstdc++)
# so a newer libstdc++ isn't shadowed by an older system one.
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})

add_subdirectory(src/cafnusyst)

install(EXPORT cafnusyst-targets
        NAMESPACE cafnusyst::
        DESTINATION lib/cmake/cafnusyst)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/setup.cafnusyst.sh.in 
  ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh @ONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh DESTINATION bin)
