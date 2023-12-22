cmake_minimum_required(VERSION 3.14)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
    message(STATUS "No build type selected, default to ${CMAKE_BUILD_TYPE}")
endif()

project(AMGCL-C C CXX)

set(PACKAGE_VERSION   0.1.0)

option(AMGCL-C_BUILD_EXAMPLES OFF)

find_package(OpenMP REQUIRED)

find_package(Boost COMPONENTS serialization)

include(ExternalProject)

# This  just downloads and unpacks amgcl to build/amgcl.
# No compilation here as we access it as a header only library.
set(AMGCL_DIR "${CMAKE_CURRENT_BINARY_DIR}/amgcl")

ExternalProject_Add(AMGCL
  URL https://github.com/ddemidov/amgcl/archive/refs/tags/1.4.4.tar.gz
  PREFIX "${AMGCL_DIR}" #  Root directory for the external project. 
  TMP_DIR "tmp"
  STAMP_DIR "stamps" # Directory in which to store the timestamps of each step
  LOG_DIR "logs" # Directory in which to store the logs of each step.
  DOWNLOAD_DIR "downloads" # Directory in which to store downloaded files before unpacking them
  SOURCE_DIR "${AMGCL_DIR}" # Source directory into which downloaded contents will be unpacked,
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" #  build directory location
  INSTALL_DIR "${AMGCL_DIR}" # Source directory into which downloaded contents will be unpacked,
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

include_directories(include)
add_subdirectory(include)
add_subdirectory(src)

if (AMGCL_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

install(TARGETS amgcl-c
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
  )

install(FILES include/amgcl-c.h
  DESTINATION include/amgcl-c
  )