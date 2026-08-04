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

set(ANAOBJ_LIBS
    duneanaobj::StandardRecordProxy
    duneanaobj::StandardRecordFlat
    duneanaobj::StandardRecord
    duneanaobj::StandardRecord_dict
)

find_package(GENIE3 REQUIRED)
if(NOT TARGET GENIE3::All)
  message(FATAL_ERROR "Expected find_package(GENIE3 REQUIRED) call to set up target GENIE3::All.")
endif()

## Check if we have fhiclcpp from ups
find_package(fhiclcpp QUIET)
if(NOT fhiclcpp_FOUND)
  find_package(fhiclcppstandalone)
endif()

###### Compiler set up
add_library(cafnusyst_dependencies INTERFACE)
target_link_libraries(cafnusyst_dependencies INTERFACE fhiclcpp::fhiclcpp)
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

add_subdirectory(src/cafnusyst)

add_subdirectory(fcls)

install(EXPORT cafnusyst-targets
        NAMESPACE cafnusyst::
        DESTINATION lib/cmake/cafnusyst)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/setup.cafnusyst.sh.in 
  ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh @ONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh DESTINATION bin)
