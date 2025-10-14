cmake_minimum_required(VERSION 3.16)

# Project metadata
project(spalt VERSION 0.1.0 LANGUAGES CXX)
message(STATUS "Current source directory: ${CMAKE_CURRENT_SOURCE_DIR}")

# Allow user to choose standard >= 17; default to 20 if available, else 17
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(SPALT_BUILD_TESTS "Build tests" ON)

# Warnings and sanitizers (optional, can be toggled by user)
if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Eigen3 CONFIG)
if(NOT Eigen3_FOUND)
  message(STATUS "System Eigen not found. Download Eigen 3.4.0.")
  include(FetchContent)
  FetchContent_Populate(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  )
  add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${eigen3_SOURCE_DIR}"
  )
endif()
message(STATUS "Eigen3 found: ${Eigen3_FOUND}")


# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find Open3D (after Eigen3)
find_package(Open3D REQUIRED)

# Find RDKit
# First try to find RDKit via CMake config (conda/system installs)
find_package(RDKit QUIET)

if(RDKit_FOUND)
  message(STATUS "RDKit found via find_package")
else()
  message(STATUS "RDKit not found via find_package, trying manual search...")

  # Manual search with improved path detection
  set(RDKIT_ROOT "" CACHE PATH "Path to RDKit installation")

  # Find RDKit include directory
  find_path(RDKIT_INCLUDE_DIR
    NAMES "rdkit/GraphMol/GraphMol.h" "GraphMol/GraphMol.h"
    PATHS
      "${RDKIT_ROOT}/include"
      "${RDKIT_ROOT}/Code"
      "$ENV{CONDA_PREFIX}/include"
      "$ENV{CONDA_PREFIX}/include/rdkit"
      "$ENV{HOME}/.local/include"
      "$ENV{HOME}/.local/include/rdkit"
      "$ENV{HOME}/.local/lib/python*/site-packages"
      "/usr/include"
      "/usr/include/rdkit"
      "/usr/local/include"
      "/usr/local/include/rdkit"
      "/opt/rdkit/include"
      "/opt/rdkit/include/rdkit"
  )

  # List of required RDKit libraries
  set(RDKIT_LIBRARIES_NEEDED
    RDKitGraphMol
    RDKitSmilesParse
    RDKitFileParsers
    RDKitRDGeneral
    RDKitSubstructMatch
    RDKitDistGeometry
    RDKitDistGeomHelpers
    RDKitPartialCharges
    RDKitForceFieldHelpers
    RDKitForceField
    RDKitChemTransforms
    RDKitDescriptors
    RDKitMolAlign
  )

  # Find each library
  foreach(lib ${RDKIT_LIBRARIES_NEEDED})
    find_library(${lib}_LIBRARY
      NAMES ${lib}
      PATHS
        "${RDKIT_ROOT}/lib"
        "${RDKIT_ROOT}/build/lib"
        "$ENV{CONDA_PREFIX}/lib"
        "$ENV{HOME}/.local/lib"
        "$ENV{HOME}/.local/lib/python*/site-packages/rdkit.libs"
        "/usr/lib"
        "/usr/local/lib"
        "/opt/rdkit/lib"
    )
    if(${lib}_LIBRARY)
      list(APPEND RDKIT_LIBRARIES ${${lib}_LIBRARY})
    else()
      message(WARNING "Could not find ${lib} library")
    endif()
  endforeach()

  # Check if we found everything
  if(NOT RDKIT_INCLUDE_DIR)
    message(FATAL_ERROR "RDKit include directory not found. Please install RDKit with development headers.\n"
      "Note: pip install rdkit does NOT include C++ headers - use conda or apt instead.\n"
      "Common installation methods:\n"
      "  conda: conda install -c conda-forge rdkit\n"
      "  apt: sudo apt install librdkit-dev rdkit-data\n"
      "  Manual: cmake -DRDKIT_ROOT=/path/to/rdkit")
  endif()

  # Check for essential libraries
  set(RDKIT_ESSENTIAL_LIBS RDKitGraphMol RDKitSmilesParse RDKitFileParsers RDKitRDGeneral)
  foreach(lib ${RDKIT_ESSENTIAL_LIBS})
    if(NOT ${lib}_LIBRARY)
      message(FATAL_ERROR "Essential RDKit library ${lib} not found.\n"
        "Found include dir: ${RDKIT_INCLUDE_DIR}\n"
        "Found libraries: ${RDKIT_LIBRARIES}")
    endif()
  endforeach()

  message(STATUS "Found RDKit include: ${RDKIT_INCLUDE_DIR}")
  message(STATUS "Found RDKit libraries: ${RDKIT_LIBRARIES}")

  # Create RDKit interface library for easier linking
  add_library(RDKit INTERFACE IMPORTED GLOBAL)
  set_target_properties(RDKit PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${RDKIT_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${RDKIT_LIBRARIES}"
  )
endif()

# Create spalt library
add_library(spalt_lib
  src/surface.cpp
  src/utils.cpp
  src/fps.cpp
  src/molecule.cpp
  src/props/esp.cpp
  src/props/hydrophobicity.cpp
  src/props/hbond.cpp
  src/props/rdkit_utils.cpp
  src/props/property_params.cpp
  src/props/pharmacophore.cpp
  src/conformer_generator.cpp
  src/surface_fibonacci.cpp
)

# Set library properties
set_target_properties(spalt_lib PROPERTIES
  OUTPUT_NAME spalt
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Link dependencies for library
target_link_libraries(spalt_lib PUBLIC
    Eigen3::Eigen
    Open3D::Open3D
    RDKit
)

target_compile_features(spalt_lib PUBLIC cxx_std_${CMAKE_CXX_STANDARD})
target_compile_definitions(spalt_lib PUBLIC SPALT_HAS_RDKIT=1)
target_include_directories(spalt_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Create executable
add_executable(spalt src/cli.cpp)
target_link_libraries(spalt PRIVATE spalt_lib)
target_compile_features(spalt PRIVATE cxx_std_${CMAKE_CXX_STANDARD})

# Install rules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install library
install(TARGETS spalt_lib
  EXPORT spalt-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install executable
install(TARGETS spalt
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install headers
install(DIRECTORY src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# Generate and install CMake config files
set(SPALT_CMAKE_CONFIG_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/spalt")

# Create spalt-config.cmake
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/spalt-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spalt-config.cmake"
  INSTALL_DESTINATION ${SPALT_CMAKE_CONFIG_DIR}
  PATH_VARS
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_LIBDIR
)

# Create spalt-config-version.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/spalt-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install CMake config files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/spalt-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/spalt-config-version.cmake"
  DESTINATION ${SPALT_CMAKE_CONFIG_DIR}
)

# Install targets file
install(EXPORT spalt-targets
  FILE spalt-targets.cmake
  NAMESPACE spalt::
  DESTINATION ${SPALT_CMAKE_CONFIG_DIR}
)

if(SPALT_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()