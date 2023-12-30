cmake_minimum_required(VERSION 3.14)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
  message(STATUS "No build type selected, default to ${CMAKE_BUILD_TYPE}")
endif()


# Set  CMAKE_EXPORT_PACKAGE_REGISTRY  for exporting via user registry.
# See https://cmake.org/cmake/help/latest/command/export.html#exporting-packages

project(AMGCL_C C CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(ExternalProject)
include(CTest)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(AMGCL_VERSION 1.4.4)
set(AMGCL_C_MAJOR_VERSION 0)
set(AMGCL_C_MINOR_VERSION 1)
set(AMGCL_C_PATCH_VERSION 0)
set(AMGCL_C_VERSION    ${AMGCL_C_MAJOR_VERSION}.${AMGCL_C_MINOR_VERSION}.${AMGCL_C_PATCH_VERSION})


find_package(Boost COMPONENTS serialization REQUIRED)

# I would very much like to REQUIRE this, but on apple
# there are problems I wasn't able to solve yet
# Any help appreciated
find_package(OpenMP)
if (OpenMP_FOUND)
  message(STATUS "Found OpenMP")
else()
  message(STATUS "Sorry, no OpenMP support")
endif()

# This just downloads and unpacks amgcl to build/amgcl.
# No compilation here as we access it as a header only library.
set(AMGCL_DIR "${CMAKE_CURRENT_BINARY_DIR}/amgcl")

ExternalProject_Add(AMGCL
  URL https://github.com/ddemidov/amgcl/archive/refs/tags/${AMGCL_VERSION}.tar.gz
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

# AMGCL_C code
include_directories(${AMGCL_DIR})
add_subdirectory(src)
add_subdirectory(examples)

#
# Export and Install
#
# From here, see https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html
install(TARGETS amgcl_c
  EXPORT AMGCL_CTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(FILES include/amgcl_c.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amgcl_c
  )


install(EXPORT AMGCL_CTargets
  FILE AMGCL_CTargets.cmake
  NAMESPACE AMGCL_C::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL_C
  )

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL_CConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL_C
  )

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL_CConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL_CConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL_C
  )



set_property(TARGET amgcl_c PROPERTY VERSION ${AMGCL_C_VERSION})
set_property(TARGET amgcl_c PROPERTY SOVERSION ${AMGCL_C_MAJOR_VERSION})
set_property(TARGET amgcl_c PROPERTY INTERFACE_AMGCL_C_MAJOR_VERSION ${AMGCL_C_MAJOR_VERSION})
set_property(TARGET amgcl_c APPEND PROPERTY COMPATIBLE_INTERFACE_STRING AMGCL_C_MAJOR_VERSION)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL_CConfigVersion.cmake"
  VERSION "${AMGCL_C_VERSION}"
  COMPATIBILITY AnyNewerVersion
  )

export(EXPORT AMGCL_CTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/AMGCL_CTargets.cmake"
       NAMESPACE AMGCL_C::
       )



     