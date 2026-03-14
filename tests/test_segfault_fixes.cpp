#include <gtest/gtest.h>
#include "conformer_params.h"
#include "molecule.h"
#include "surface.h"

/**
 * @brief Test suite for segmentation fault fixes
 *
 * Tests edge cases and boundary conditions that previously caused segmentation faults.
 * These tests verify that bounds checking and safety measures work correctly.
 */
class SegfaultFixesTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Set up test molecule
        mol = std::make_unique<Molecule>("CCO");  // Simple ethanol molecule
    }

    std::unique_ptr<Molecule> mol;
};

TEST_F(SegfaultFixesTest, ConformerGeneratorVariableCollisionFix) {
    // Test that the variable name collision in conformer generator is fixed
    ConformerParams params;
    params.total_conformers = 20;
    params.num_clusters = 5;
    params.enable_optimization = true;

    // This should not cause an infinite loop or crash
    auto cluster_map = mol->generateAdvancedConformers(params);

    // Should have some conformers (may be fewer than requested due to clustering)
    EXPECT_GE(mol->getNumConformers(), 1);
    EXPECT_LE(mol->getNumConformers(), 5);  // Should not exceed num_clusters

    // All conformers should be valid
    for (int i = 0; i < mol->getNumConformers(); ++i) {
        EXPECT_NO_THROW({
            auto& conf = mol->get_mol()->getConformer(i);
            EXPECT_TRUE(conf.is3D());
        });
    }
}

TEST_F(SegfaultFixesTest, ConformerAccessBoundsCheck) {
    // Test bounds checking for conformer access
    ConformerParams params;
    params.total_conformers = 10;
    params.num_clusters = 3;

    auto cluster_map = mol->generateAdvancedConformers(params);

    // Verify that all conformer IDs in cluster_map are valid
    for (const auto& [cluster_id, conformer_id] : cluster_map) {
        EXPECT_GE(conformer_id, 0);
        EXPECT_LT(conformer_id, mol->get_mol()->getNumConformers());

        // Should be able to access conformer without crashing
        EXPECT_NO_THROW({
            auto& conf = mol->get_mol()->getConformer(conformer_id);
            EXPECT_TRUE(conf.is3D());
        });
    }
}

TEST_F(SegfaultFixesTest, SurfaceVertexBoundsCheck) {
    // Test vertex bounds checking in surface generation
    SurfaceParams params;
    params.num_vertices = 50;
    params.program = "msms";  // Use string, not enum

    // Ensure conformers exist first
    mol->calculate_conformers(1);
    EXPECT_GT(mol->getNumConformers(), 0);

    // Generate surface using createSurface
    int surface_id =
        mol->createSurface(0, params.num_vertices, params.radius, params.density, params.hdensity,
                           params.type, params.program, params.sample_method);
    EXPECT_GE(surface_id, 0);

    // Get surface vertices for testing
    auto vertices = mol->getSurfaceVertices(0);
    EXPECT_GT(vertices.size(), 0);
}

TEST_F(SegfaultFixesTest, SurfaceVertexNormalComputationBounds) {
    // Test bounds checking in vertex normal computation
    SurfaceParams params;
    params.num_vertices = 50;
    params.program = "fibonacci";  // Use string, not enum

    // Ensure conformers exist first
    mol->calculate_conformers(1);
    EXPECT_GT(mol->getNumConformers(), 0);

    // Generate Fibonacci surface (which computes normals)
    int surface_id =
        mol->createSurface(0, params.num_vertices, params.radius, params.density, params.hdensity,
                           params.type, params.program, params.sample_method);
    EXPECT_GE(surface_id, 0);

    // Get surface vertices for testing
    auto vertices = mol->getSurfaceVertices(0);
    EXPECT_GT(vertices.size(), 0);

    // Verify that all vertices have proper normals (6 elements: x,y,z,nx,ny,nz)
    for (const auto& [index, vertex_data] : vertices) {
        const auto& coords = vertex_data.second;
        EXPECT_GE(coords.size(), 6) << "Vertex " << index << " has insufficient coordinates";

        if (coords.size() >= 6) {
            // Check that normal components are reasonable (not NaN or infinite)
            EXPECT_FALSE(std::isnan(coords[3])) << "Normal x is NaN for vertex " << index;
            EXPECT_FALSE(std::isnan(coords[4])) << "Normal y is NaN for vertex " << index;
            EXPECT_FALSE(std::isnan(coords[5])) << "Normal z is NaN for vertex " << index;

            EXPECT_FALSE(std::isinf(coords[3])) << "Normal x is infinite for vertex " << index;
            EXPECT_FALSE(std::isinf(coords[4])) << "Normal y is infinite for vertex " << index;
            EXPECT_FALSE(std::isinf(coords[5])) << "Normal z is infinite for vertex " << index;
        }
    }
}

TEST_F(SegfaultFixesTest, EdgeCaseEmptyConformerList) {
    // Test handling of edge case where no conformers are generated
    ConformerParams params;
    params.total_conformers = 0;  // Force no conformers
    params.num_clusters = 3;

    // This should not crash
    auto cluster_map = mol->generateAdvancedConformers(params);

    // Should handle gracefully
    EXPECT_TRUE(cluster_map.empty());
    EXPECT_EQ(mol->getNumConformers(), 0);
}

TEST_F(SegfaultFixesTest, EdgeCaseMoreClustersThanConformers) {
    // Test edge case where num_clusters > total_conformers
    ConformerParams params;
    params.total_conformers = 2;
    params.num_clusters = 5;  // More clusters than conformers

    // This should not crash
    auto cluster_map = mol->generateAdvancedConformers(params);

    // Should return all available conformers
    EXPECT_LE(cluster_map.size(), 2);
    EXPECT_LE(mol->getNumConformers(), 2);
}

TEST_F(SegfaultFixesTest, LargeConformerGeneration) {
    // Test with larger number of conformers to stress test the fixes
    ConformerParams params;
    params.total_conformers = 100;
    params.num_clusters = 10;
    params.enable_optimization = true;

    // This should not cause infinite loops or crashes
    EXPECT_NO_THROW({
        auto cluster_map = mol->generateAdvancedConformers(params);

        // Should have reasonable number of conformers
        EXPECT_GE(mol->getNumConformers(), 1);
        EXPECT_LE(mol->getNumConformers(), 10);

        // All conformers should be valid
        for (int i = 0; i < mol->getNumConformers(); ++i) {
            auto& conf = mol->get_mol()->getConformer(i);
            EXPECT_TRUE(conf.is3D());
        }
    });
}

TEST_F(SegfaultFixesTest, SurfaceSubsamplingBounds) {
    // Test surface subsampling with various edge cases
    SurfaceParams params;
    params.num_vertices = 1000;  // Request more vertices than available
    params.program = "msms";     // Use string, not enum

    // Ensure conformers exist first
    mol->calculate_conformers(1);
    EXPECT_GT(mol->getNumConformers(), 0);

    // Generate surface using createSurface
    int surface_id =
        mol->createSurface(0, params.num_vertices, params.radius, params.density, params.hdensity,
                           params.type, params.program, params.sample_method);
    EXPECT_GE(surface_id, 0);

    // Get surface vertices for testing
    auto vertices = mol->getSurfaceVertices(0);
    EXPECT_GT(vertices.size(), 0);

    // Test that we can access vertices without bounds errors
    EXPECT_NO_THROW({
        EXPECT_GT(vertices.size(), 0);
        // Check a few vertices have proper structure
        for (const auto& [index, vertex_data] : vertices) {
            if (index >= 5)
                break;  // Just check first few
            EXPECT_GE(vertex_data.second.size(), 6);
        }
    });
}
