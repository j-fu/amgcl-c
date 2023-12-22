cmake_minimum_required(VERSION 3.14)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
    message(STATUS "No build type selected, default to ${CMAKE_BUILD_TYPE}")
endif()

# Set  CMAKE_EXPORT_PACKAGE_REGISTRY  for exporting via user registry.
# See https://cmake.org/cmake/help/latest/command/export.html#exporting-packages

project(AMGCL-C C CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(AMGCL-C_MAJOR_VERSION 0)
set(AMGCL-C_MINOR_VERSION 1)
set(AMGCL-C_PATCH_VERSION 0)
set(AMGCL-C_VERSION    ${AMGCL-C_MAJOR_VERSION}.${AMGCL-C_MINOR_VERSION}.${AMGCL-C_PATCH_VERSION})


option(AMGCL-C_BUILD_EXAMPLES OFF)

find_package(OpenMP REQUIRED)

find_package(Boost COMPONENTS serialization)

include(ExternalProject)

# This  just downloads and unpacks amgcl to build/amgcl.
# No compilation here as we access it as a header only library.
# when we compile amgcl-c
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

include_directories(${AMGCL_DIR})

add_subdirectory(src)

if (AMGCL-C_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# From here, see https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html
install(TARGETS amgcl-c
  EXPORT AMGCL-CTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(FILES include/amgcl-c.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amgcl-c
  )


install(EXPORT AMGCL-CTargets
  FILE AMGCL-CTargets.cmake
  NAMESPACE AMGCL-C::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL-C
  )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL-CConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL-C
  )

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL-CConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL-CConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMGCL-C
  )



set_property(TARGET amgcl-c PROPERTY VERSION ${AMGCL-C_VERSION})
set_property(TARGET amgcl-c PROPERTY SOVERSION ${AMGCL-C_MAJOR_VERSION})
set_property(TARGET amgcl-c PROPERTY
  INTERFACE_AMGCL-C_MAJOR_VERSION ${AMGCL-C_MAJOR_VERSION})
set_property(TARGET amgcl-c APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING AMGCL-C_MAJOR_VERSION
)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/AMGCL-CConfigVersion.cmake"
  VERSION "${AMGCL-C_VERSION}"
  COMPATIBILITY AnyNewerVersion
  )

export(EXPORT AMGCL-CTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/AMGCL-CTargets.cmake"
       NAMESPACE AMGCL-C::
       )

     