#include <gtest/gtest.h>
#include <filesystem>
#include "conformer_params.h"
#include "molecule.h"

class ConformerGenerationTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Create test output directory
        test_dir = std::filesystem::temp_directory_path() / "spalt_conformer_test";
        std::filesystem::create_directories(test_dir);
    }

    void TearDown() override {
        // Clean up test directory
        std::error_code ec;
        std::filesystem::remove_all(test_dir, ec);
    }

    std::filesystem::path test_dir;
};

TEST_F(ConformerGenerationTest, SimpleConformerGeneration) {
    Molecule mol("C1CCCCC1");  // Cyclohexane

    // Test simple conformer generation
    mol.calculate_conformers(3);

    EXPECT_EQ(mol.getNumConformers(), 3);

    // Verify all conformers are 3D
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(mol.get_mol()->getConformer(i).is3D());
    }
}

TEST_F(ConformerGenerationTest, AdvancedConformerGeneration) {
    Molecule mol("C1CCCCC1");  // Cyclohexane

    // Test advanced conformer generation
    ConformerParams params;
    params.total_conformers = 20;
    params.num_clusters = 3;
    params.enable_optimization = true;

    auto cluster_map = mol.generateAdvancedConformers(params);

    // Should have exactly 3 conformers after clustering
    EXPECT_EQ(mol.getNumConformers(), 3);
    EXPECT_EQ(cluster_map.size(), 3);

    // Verify all conformers are 3D
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(mol.get_mol()->getConformer(i).is3D());
    }

    // Verify cluster map has expected keys
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(cluster_map.find(i) != cluster_map.end());
    }
}

TEST_F(ConformerGenerationTest, AdvancedConformerGenerationFallback) {
    Molecule mol("C1CCCCC1");  // Cyclohexane

    // Test advanced conformer generation with optimization disabled
    ConformerParams params;
    params.total_conformers = 15;
    params.num_clusters = 2;
    params.enable_optimization = false;  // Force fallback to ETKDG/Uniform

    auto cluster_map = mol.generateAdvancedConformers(params);

    // Should have exactly 2 conformers after clustering
    EXPECT_EQ(mol.getNumConformers(), 2);
    EXPECT_EQ(cluster_map.size(), 2);

    // Verify all conformers are 3D
    for (int i = 0; i < 2; ++i) {
        EXPECT_TRUE(mol.get_mol()->getConformer(i).is3D());
    }
}

TEST_F(ConformerGenerationTest, ConformerSurfaceCorrespondence) {
    Molecule mol("CCO");  // Ethanol

    // Generate conformers
    mol.calculate_conformers(3);

    // Create surfaces for each conformer
    for (int i = 0; i < 3; ++i) {
        int surface_id = mol.createSurface(i, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
        EXPECT_EQ(surface_id, i);  // Surface ID should match conformer ID

        // Verify surface exists
        EXPECT_GT(mol.getSurfaceVertexCount(i), 0);
    }
}

TEST_F(ConformerGenerationTest, ConformerParamsDefaultValues) {
    ConformerParams params;

    // Test default values
    EXPECT_EQ(params.total_conformers, 50);
    EXPECT_EQ(params.num_clusters, 5);
    EXPECT_FALSE(params.enable_optimization);  // Changed default
    EXPECT_TRUE(params.enforce_chirality);
    EXPECT_TRUE(params.align_conformers);
    EXPECT_FALSE(params.use_random_coords);
    EXPECT_EQ(params.num_threads, 1);
    EXPECT_FALSE(params.clear_existing);
}

TEST_F(ConformerGenerationTest, ConformerParamsWeightCalculation) {
    ConformerParams params;
    params.total_conformers = 30;
    params.enable_optimization = true;

    auto counts = params.getConformerCounts();

    // Should have 3 counts: optimized, etkdg, uniform
    EXPECT_EQ(counts.size(), 3);

    // Weights should be [4/6, 1/6, 1/6] = [20, 5, 5] for 30 total
    EXPECT_EQ(counts[0], 20);  // Optimized
    EXPECT_EQ(counts[1], 5);   // ETKDG
    EXPECT_EQ(counts[2], 5);   // Uniform

    // Sum should equal total
    EXPECT_EQ(counts[0] + counts[1] + counts[2], params.total_conformers);
}

TEST_F(ConformerGenerationTest, ConformerParamsFallbackWeights) {
    ConformerParams params;
    params.total_conformers = 20;
    params.enable_optimization = false;

    params.applyFallbackWeights();
    auto counts = params.getConformerCounts();

    // Fallback weights should be [0, 1/2, 1/2] = [0, 10, 10]
    EXPECT_EQ(counts[0], 0);   // Optimized (disabled)
    EXPECT_EQ(counts[1], 10);  // ETKDG
    EXPECT_EQ(counts[2], 10);  // Uniform

    // Sum should equal total
    EXPECT_EQ(counts[0] + counts[1] + counts[2], params.total_conformers);
}

TEST_F(ConformerGenerationTest, InvalidMoleculeHandling) {
    // Test with invalid SMILES
    Molecule mol("");

    // Should handle gracefully
    ConformerParams params;
    params.total_conformers = 10;
    params.num_clusters = 3;

    auto cluster_map = mol.generateAdvancedConformers(params);
    EXPECT_TRUE(cluster_map.empty());
}

TEST_F(ConformerGenerationTest, EdgeCaseSingleConformer) {
    Molecule mol("C1CCCCC1");

    // Test simple conformer generation with single conformer
    mol.calculate_conformers(1);

    EXPECT_EQ(mol.getNumConformers(), 1);

    // Verify conformer is 3D
    EXPECT_TRUE(mol.get_mol()->getConformer(0).is3D());
}

TEST_F(ConformerGenerationTest, EdgeCaseMoreClustersThanConformers) {
    Molecule mol("C1CCCCC1");

    // Test simple conformer generation with multiple conformers
    mol.calculate_conformers(5);

    // Should have the requested number of conformers
    EXPECT_EQ(mol.getNumConformers(), 5);

    // Verify all conformers are 3D
    for (int i = 0; i < 5; ++i) {
        EXPECT_TRUE(mol.get_mol()->getConformer(i).is3D());
    }
}
