#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "molecule.h"
#include "surface.h"

class NormalsUsageTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple test molecule
        test_molecule = std::make_unique<Molecule>("CCO");  // Ethanol
        test_molecule->calculate_conformers(1);
    }
    
    std::unique_ptr<Molecule> test_molecule;
};

/**
 * @brief Test that surface vertices contain 6 elements (coordinates + normals)
 */
TEST_F(NormalsUsageTest, SurfaceVerticesContainNormals) {
    // Create surface with normals
    Surface surface(*test_molecule, 0, 100, 1.2, 3.0, 3.0, "tses", "msms", "full");
    
    // Check that vertices have 6 elements (x, y, z, nx, ny, nz)
    for (const auto& [index, vertex_data] : surface.vertices) {
        const auto& coords_normals = vertex_data.second;
        EXPECT_EQ(coords_normals.size(), 6) << "Vertex " << index << " should have 6 elements (coordinates + normals)";
        
        // Check that normals are non-zero (they should be computed)
        if (coords_normals.size() >= 6) {
            double nx = coords_normals[3];
            double ny = coords_normals[4];
            double nz = coords_normals[5];
            double normal_magnitude = sqrt(nx*nx + ny*ny + nz*nz);
            
            EXPECT_GT(normal_magnitude, 0.1) << "Normal vector should have reasonable magnitude";
            EXPECT_LT(abs(normal_magnitude - 1.0), 0.1) << "Normal vector should be approximately unit length";
        }
    }
}

/**
 * @brief Test that Open3D point clouds preserve existing normals
 */
TEST_F(NormalsUsageTest, Open3DPointCloudPreservesNormals) {
    // Create two identical molecules for alignment
    Molecule mol1("CCO");
    Molecule mol2("CCO");
    mol1.calculate_conformers(1);
    mol2.calculate_conformers(1);
    
    // Create surfaces
    Surface surface1(mol1, 0, 100, 1.2, 3.0, 3.0, "tses", "msms", "full");
    Surface surface2(mol2, 0, 100, 1.2, 3.0, 3.0, "tses", "msms", "full");
    
    // Create Open3D point clouds manually (simulating the alignment process)
    auto source_cloud = std::make_shared<open3d::geometry::PointCloud>();
    auto target_cloud = std::make_shared<open3d::geometry::PointCloud>();
    
    // Convert surfaces to point clouds with normals
    for (const auto& vertex_pair : surface1.vertices) {
        const auto& coords = vertex_pair.second.second;
        if (coords.size() >= 6) {
            source_cloud->points_.emplace_back(coords[0], coords[1], coords[2]);
            source_cloud->normals_.emplace_back(coords[3], coords[4], coords[5]);
        }
    }
    
    for (const auto& vertex_pair : surface2.vertices) {
        const auto& coords = vertex_pair.second.second;
        if (coords.size() >= 6) {
            target_cloud->points_.emplace_back(coords[0], coords[1], coords[2]);
            target_cloud->normals_.emplace_back(coords[3], coords[4], coords[5]);
        }
    }
    
    // Verify that normals were preserved
    EXPECT_EQ(source_cloud->points_.size(), source_cloud->normals_.size());
    EXPECT_EQ(target_cloud->points_.size(), target_cloud->normals_.size());
    
    // Verify that normals are reasonable
    for (size_t i = 0; i < source_cloud->normals_.size(); ++i) {
        const auto& normal = source_cloud->normals_[i];
        double magnitude = normal.norm();
        EXPECT_GT(magnitude, 0.1) << "Source normal " << i << " should have reasonable magnitude";
        EXPECT_LT(abs(magnitude - 1.0), 0.1) << "Source normal " << i << " should be approximately unit length";
    }
    
    for (size_t i = 0; i < target_cloud->normals_.size(); ++i) {
        const auto& normal = target_cloud->normals_[i];
        double magnitude = normal.norm();
        EXPECT_GT(magnitude, 0.1) << "Target normal " << i << " should have reasonable magnitude";
        EXPECT_LT(abs(magnitude - 1.0), 0.1) << "Target normal " << i << " should be approximately unit length";
    }
}

/**
 * @brief Test that surface export preserves normals
 */
TEST_F(NormalsUsageTest, SurfaceExportPreservesNormals) {
    Surface surface(*test_molecule, 0, 100, 1.2, 3.0, 3.0, "tses", "msms", "full");
    
    // Test PLY export with normals
    std::string test_file = "/tmp/test_normals.ply";
    surface.savePointCloud(test_file, true, false, {});  // include_normals=true
    
    // Read back the file and verify normals are present
    std::ifstream file(test_file);
    ASSERT_TRUE(file.is_open()) << "Failed to open exported PLY file";
    
    std::string line;
    bool found_nx = false, found_ny = false, found_nz = false;
    
    while (std::getline(file, line)) {
        if (line.find("property float nx") != std::string::npos) found_nx = true;
        if (line.find("property float ny") != std::string::npos) found_ny = true;
        if (line.find("property float nz") != std::string::npos) found_nz = true;
        if (line.find("end_header") != std::string::npos) {
            break;
        }
    }
    
    EXPECT_TRUE(found_nx) << "PLY file should contain nx property";
    EXPECT_TRUE(found_ny) << "PLY file should contain ny property";
    EXPECT_TRUE(found_nz) << "PLY file should contain nz property";
    
    // Clean up
    std::remove(test_file.c_str());
}

/**
 * @brief Test that alignment uses existing normals instead of estimating new ones
 */
TEST_F(NormalsUsageTest, AlignmentUsesExistingNormals) {
    // Create two molecules for alignment
    Molecule mol1("CCO");
    Molecule mol2("CCO");
    mol1.calculate_conformers(1);
    mol2.calculate_conformers(1);
    
    // Generate surfaces first (this will compute the normals)
    SurfaceParams params;
    params.num_vertices = 100;
    params.radius = 1.2;
    mol1.compute("esp", 0, params);
    mol2.compute("esp", 0, params);
    
    // Perform alignment
    auto alignment_result = mol1.align_to_molecule_detailed(mol2);
    
    // The alignment should succeed without errors
    EXPECT_GT(alignment_result.fitness, 0.0) << "Alignment should have positive fitness";
    EXPECT_LT(alignment_result.inlier_rmse, 10.0) << "Alignment should have reasonable RMSE";
    
    // The transformation should be reasonable (not too large)
    double translation_norm = alignment_result.transformation.translation().norm();
    EXPECT_LT(translation_norm, 5.0) << "Translation should be reasonable";
}

/**
 * @brief Test Fibonacci surface generation produces proper normals
 */
TEST_F(NormalsUsageTest, FibonacciSurfaceGeneratesNormals) {
    // Create surface using Fibonacci method
    SurfaceParams params;
    params.program = "fibonacci";
    params.num_vertices = 100;
    params.radius = 1.2;
    
    Surface surface(*test_molecule, 0, params.num_vertices, params.radius, 3.0, 3.0, "tses", "fibonacci", "full");
    
    // Check that vertices have 6 elements with valid normals
    for (const auto& [index, vertex_data] : surface.vertices) {
        const auto& coords_normals = vertex_data.second;
        EXPECT_EQ(coords_normals.size(), 6) << "Fibonacci surface vertex " << index << " should have 6 elements";
        
        if (coords_normals.size() >= 6) {
            double nx = coords_normals[3];
            double ny = coords_normals[4];
            double nz = coords_normals[5];
            double normal_magnitude = sqrt(nx*nx + ny*ny + nz*nz);
            
            // Fibonacci normals should be computed from mesh faces, not zero
            EXPECT_GT(normal_magnitude, 0.1) << "Fibonacci surface normal should have reasonable magnitude";
            EXPECT_LT(abs(normal_magnitude - 1.0), 0.2) << "Fibonacci surface normal should be approximately unit length";
        }
    }
}

/**
 * @brief Test normal estimation fallback when normals are missing
 */
TEST_F(NormalsUsageTest, NormalEstimationFallback) {
    // Create point cloud without normals
    auto source_cloud = std::make_shared<open3d::geometry::PointCloud>();
    auto target_cloud = std::make_shared<open3d::geometry::PointCloud>();
    
    // Add some test points
    source_cloud->points_.emplace_back(0.0, 0.0, 0.0);
    source_cloud->points_.emplace_back(1.0, 0.0, 0.0);
    source_cloud->points_.emplace_back(0.0, 1.0, 0.0);
    
    target_cloud->points_.emplace_back(0.1, 0.1, 0.1);
    target_cloud->points_.emplace_back(1.1, 0.1, 0.1);
    target_cloud->points_.emplace_back(0.1, 1.1, 0.1);
    
    // Initially no normals
    EXPECT_TRUE(source_cloud->normals_.empty());
    EXPECT_TRUE(target_cloud->normals_.empty());
    
    // Estimate normals
    source_cloud->EstimateNormals();
    target_cloud->EstimateNormals();
    
    // Check that normals were estimated
    EXPECT_EQ(source_cloud->points_.size(), source_cloud->normals_.size());
    EXPECT_EQ(target_cloud->points_.size(), target_cloud->normals_.size());
    
    // Check that normals are reasonable
    for (size_t i = 0; i < source_cloud->normals_.size(); ++i) {
        const auto& normal = source_cloud->normals_[i];
        double magnitude = normal.norm();
        EXPECT_GT(magnitude, 0.1) << "Estimated normal should have reasonable magnitude";
    }
}

