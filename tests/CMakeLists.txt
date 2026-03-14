include(FetchContent)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
else()
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  )
endif()

# Recommended by GoogleTest on Windows; harmless elsewhere
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(spalt_tests
  test_surface.cpp
  test_utils.cpp
  test_molecule.cpp
  test_esp.cpp
  test_hbond_features.cpp
  test_hydrophobicity.cpp
  test_conformer_generation.cpp
  test_property_system.cpp
  test_alignment.cpp
  test_fibonacci_surface.cpp
  test_cli_integration.cpp
  test_normals_usage.cpp
  test_segfault_fixes.cpp
  test_tsv_parsing.cpp
  test_pharmacophore.cpp
)

# Give the tests access to headers in src/
target_include_directories(spalt_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Compile Surface implementation into the test binary
 target_sources(spalt_tests PRIVATE
   ${CMAKE_SOURCE_DIR}/src/surface.cpp
   ${CMAKE_SOURCE_DIR}/src/utils.cpp
   ${CMAKE_SOURCE_DIR}/src/fps.cpp
   ${CMAKE_SOURCE_DIR}/src/molecule.cpp
   ${CMAKE_SOURCE_DIR}/src/props/esp.cpp
   ${CMAKE_SOURCE_DIR}/src/props/hydrophobicity.cpp
   ${CMAKE_SOURCE_DIR}/src/props/hbond.cpp
   ${CMAKE_SOURCE_DIR}/src/props/pharmacophore.cpp
   ${CMAKE_SOURCE_DIR}/src/props/rdkit_utils.cpp
   ${CMAKE_SOURCE_DIR}/src/props/property_params.cpp
   ${CMAKE_SOURCE_DIR}/src/conformer_generator.cpp
   ${CMAKE_SOURCE_DIR}/src/surface_fibonacci.cpp
   ${CMAKE_SOURCE_DIR}/src/cli.cpp
 )

target_link_libraries(spalt_tests PRIVATE
    GTest::gtest_main
    Eigen3::Eigen
    Open3D::Open3D
    RDKit
)
target_compile_features(spalt_tests PRIVATE cxx_std_${CMAKE_CXX_STANDARD})
target_compile_definitions(spalt_tests PRIVATE SPALT_TESTS=1)

include(GoogleTest)
gtest_discover_tests(spalt_tests)
