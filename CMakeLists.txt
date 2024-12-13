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

# Define an option to specify the target flag
set(CAF_TYPE "" CACHE STRING "Specify the target flag: 'sbn' or 'dune'")

# Ensure the user provides a valid CAF_TYPE
if(NOT (CAF_TYPE STREQUAL "sbn" OR CAF_TYPE STREQUAL "dune"))
    message(FATAL_ERROR "CAF_TYPE must be either 'sbn' or 'dune'")
endif()

# Conditionally find the required package based on the CAF_TYPE
if(CAF_TYPE STREQUAL "sbn")

  find_package(sbnanaobj 09.20.06.03 REQUIRED)

  message(STATUS "Building against sbnanaobj")

  set(ANAOBJ_LIBS
    sbnanaobj::StandardRecord
    sbnanaobj::StandardRecordProxy
    sbnanaobj::StandardRecordFlat
  )

  add_compile_definitions(USE_SBNCAF)

elseif(CAF_TYPE STREQUAL "dune")

  find_package(duneanaobj 03.05.00 REQUIRED)

  message(STATUS "Building against duneanaobj")

  set(ANAOBJ_LIBS
    duneanaobj::StandardRecordProxy
    duneanaobj::StandardRecordFlat
    duneanaobj::StandardRecord
    duneanaobj::StandardRecord_dict
  )

  add_compile_definitions(USE_DUNECAF)

endif()



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

add_subdirectory(src/cafnusyst)

add_subdirectory(fcls)

install(EXPORT cafnusyst-targets
        NAMESPACE cafnusyst::
        DESTINATION lib/cmake/cafnusyst)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/setup.cafnusyst.sh.in 
  ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh @ONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/setup.cafnusyst.sh DESTINATION bin)
