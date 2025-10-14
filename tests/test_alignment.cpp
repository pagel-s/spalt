#include <gtest/gtest.h>
#include <filesystem>
#include "molecule.h"

class AlignmentTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Create test molecules
        ref_mol = std::make_unique<Molecule>("C1CCCCC1");  // Cyclohexane (reference)
        input_mol = std::make_unique<Molecule>("CCO");     // Ethanol (input)

        // Generate conformers
        ref_mol->calculate_conformers(1);
        input_mol->calculate_conformers(1);

        // Create surfaces
        ref_mol->createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
        input_mol->createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");

        // Create test output directory
        test_dir = std::filesystem::temp_directory_path() / "spalt_alignment_test";
        std::filesystem::create_directories(test_dir);
    }

    void TearDown() override {
        // Clean up test directory
        std::error_code ec;
        std::filesystem::remove_all(test_dir, ec);
    }

    std::unique_ptr<Molecule> ref_mol;
    std::unique_ptr<Molecule> input_mol;
    std::filesystem::path test_dir;
};

TEST_F(AlignmentTest, BasicAlignment) {
    // Test basic alignment without properties
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Verify alignment result
    EXPECT_GE(result.fitness, 0.0);
    EXPECT_LE(result.fitness, 1.0);
    EXPECT_GT(result.inlier_rmse, 0.0);
    // Note: colored ICP may not always be used depending on properties available
}

TEST_F(AlignmentTest, SinglePropertyAlignment) {
    // Compute ESP property
    input_mol->compute("esp", 0);
    ref_mol->compute("esp", 0);

    // Test alignment with single property
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Verify alignment result
    EXPECT_GT(result.fitness, 0.0);
    EXPECT_LT(result.fitness, 1.0);
    EXPECT_GT(result.inlier_rmse, 0.0);
    EXPECT_TRUE(result.used_colored_icp);
}

TEST_F(AlignmentTest, MultiPropertyAlignment) {
    // Compute multiple properties
    std::vector<std::string> properties = {"esp", "hydrophobicity", "hbond"};

    for (const auto& property : properties) {
        input_mol->compute(property, 0);
        ref_mol->compute(property, 0);
    }

    // Test alignment with multiple properties
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Verify alignment result
    EXPECT_GT(result.fitness, 0.0);
    EXPECT_LT(result.fitness, 1.0);
    EXPECT_GT(result.inlier_rmse, 0.0);
    EXPECT_TRUE(result.used_colored_icp);
}

TEST_F(AlignmentTest, MultiConformerAlignment) {
    // Generate multiple conformers for input molecule
    input_mol->calculate_conformers(3);

    // Create surfaces for all conformers
    for (int i = 0; i < 3; ++i) {
        input_mol->createSurface(i, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
    }

    // Compute properties for all conformers
    std::vector<std::string> properties = {"esp", "hydrophobicity"};

    for (int i = 0; i < 3; ++i) {
        for (const auto& property : properties) {
            input_mol->compute(property, i);
        }
    }

    // Test alignment for each conformer
    std::vector<Molecule::AlignmentResult> results;
    for (int i = 0; i < 3; ++i) {
        auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, i);
        results.push_back(result);

        // Verify alignment result
        EXPECT_GT(result.fitness, 0.0);
        EXPECT_LT(result.fitness, 1.0);
        EXPECT_GT(result.inlier_rmse, 0.0);
        EXPECT_TRUE(result.used_colored_icp);
    }

    // Results should be different for different conformers
    EXPECT_NE(results[0].fitness, results[1].fitness);
    EXPECT_NE(results[1].fitness, results[2].fitness);
}

TEST_F(AlignmentTest, ConformerSurfaceCorrespondence) {
    // Generate multiple conformers
    input_mol->calculate_conformers(3);

    // Create surfaces for all conformers
    std::vector<int> surface_ids;
    for (int i = 0; i < 3; ++i) {
        int surface_id = input_mol->createSurface(i, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
        surface_ids.push_back(surface_id);

        // Verify surface ID matches conformer ID
        EXPECT_EQ(surface_id, i);
    }

    // Compute properties and align each conformer
    std::vector<std::string> properties = {"esp"};

    for (int i = 0; i < 3; ++i) {
        // Compute properties for this conformer
        for (const auto& property : properties) {
            input_mol->compute(property, i);
        }

        // Align this conformer
        auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, i);

        // Verify alignment used correct conformer
        EXPECT_GT(result.fitness, 0.0);
        EXPECT_TRUE(result.used_colored_icp);
    }
}

TEST_F(AlignmentTest, AlignmentTransformation) {
    // Compute properties
    input_mol->compute("esp", 0);
    ref_mol->compute("esp", 0);

    // Perform alignment
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Verify transformation matrix is valid
    EXPECT_TRUE(result.transformation.matrix().allFinite());

    // Transformation should not be identity (molecules are different)
    Eigen::Matrix4d identity = Eigen::Matrix4d::Identity();
    EXPECT_FALSE(result.transformation.matrix().isApprox(identity, 1e-6));
}

TEST_F(AlignmentTest, SurfaceTransformation) {
    // Compute properties
    input_mol->compute("esp", 0);
    ref_mol->compute("esp", 0);

    // Get initial surface vertices
    auto initial_vertices = input_mol->getSurfaceVertices(0);

    // Perform alignment
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Transform surface
    EXPECT_TRUE(input_mol->transformSurface(0, result.transformation));

    // Get transformed surface vertices
    auto transformed_vertices = input_mol->getSurfaceVertices(0);

    // Vertices should be different after transformation
    EXPECT_NE(initial_vertices.size(), 0);
    EXPECT_EQ(initial_vertices.size(), transformed_vertices.size());

    // At least some vertices should be different
    bool vertices_changed = false;
    for (const auto& [atom_id, vertex_data] : initial_vertices) {
        const auto& [atom_name, initial_pos] = vertex_data;
        const auto& [_, transformed_pos] = transformed_vertices[atom_id];

        if (!initial_pos.empty() && !transformed_pos.empty()) {
            for (size_t i = 0; i < initial_pos.size(); i += 3) {
                if (i + 2 < initial_pos.size() && i + 2 < transformed_pos.size()) {
                    Eigen::Vector3d initial(initial_pos[i], initial_pos[i + 1], initial_pos[i + 2]);
                    Eigen::Vector3d transformed(transformed_pos[i], transformed_pos[i + 1],
                                                transformed_pos[i + 2]);

                    if (!initial.isApprox(transformed, 1e-6)) {
                        vertices_changed = true;
                        break;
                    }
                }
            }
        }
        if (vertices_changed)
            break;
    }

    EXPECT_TRUE(vertices_changed);
}

TEST_F(AlignmentTest, EdgeCaseSameMolecule) {
    // Test alignment of same molecule to itself
    Molecule same_mol("C1CCCCC1");
    same_mol.calculate_conformers(1);
    same_mol.createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
    same_mol.compute("esp", 0);

    auto result = same_mol.align_to_molecule_detailed(*ref_mol, 0, 0);

    // Should have good alignment (but may not be perfect due to different conformers)
    EXPECT_GT(result.fitness, 0.5);      // Good fitness for same molecule type
    EXPECT_GT(result.inlier_rmse, 0.0);  // RMSE should be positive
}

TEST_F(AlignmentTest, EdgeCaseNoProperties) {
    // Test alignment without any properties (should fall back to geometric ICP)
    auto result = input_mol->align_to_molecule_detailed(*ref_mol, 0, 0);

    // Should still work
    EXPECT_GE(result.fitness, 0.0);
    EXPECT_LE(result.fitness, 1.0);
    EXPECT_GT(result.inlier_rmse, 0.0);
    // May or may not use colored ICP depending on implementation
}

TEST_F(AlignmentTest, InvalidConformerHandling) {
    // Test with invalid conformer IDs
    EXPECT_THROW(input_mol->align_to_molecule_detailed(*ref_mol, 999, 0), std::runtime_error);
    EXPECT_THROW(input_mol->align_to_molecule_detailed(*ref_mol, 0, 999), std::runtime_error);
}
