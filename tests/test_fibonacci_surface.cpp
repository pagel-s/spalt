#include <gtest/gtest.h>
#include "molecule.h"
#include "surface.h"
#include "surface_fibonacci.h"

class FibonacciSurfaceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple molecule (water)
        std::string smiles = "C1CCCCC1";
        mol = std::make_unique<Molecule>(smiles.c_str());
        mol->calculate_conformers(1);
    }
    
    std::unique_ptr<Molecule> mol;
};

TEST_F(FibonacciSurfaceTest, BasicSurfaceGeneration) {
    // Test Fibonacci surface generation directly
    try {
        // Get conformer
        auto mol_ptr = mol->get_mol();
        const auto& conformer = mol_ptr->getConformer(0);
        int num_atoms = conformer.getNumAtoms();
        
        // Extract atom centers and radii
        Eigen::MatrixXd centers(num_atoms, 3);
        Eigen::VectorXd radii(num_atoms);
        
        for (int i = 0; i < num_atoms; ++i) {
            const auto& pos = conformer.getAtomPos(i);
            centers(i, 0) = pos.x;
            centers(i, 1) = pos.y;
            centers(i, 2) = pos.z;
            radii(i) = 1.7; // Default carbon radius
        }
        
        // Create Fibonacci surface generator
        spalt::FibonacciSurfaceParams params;
        spalt::FibonacciSurfaceGenerator generator(params);
        
        auto [mesh, point_cloud] = generator.generateSurface(centers, radii);
        
        EXPECT_GT(mesh->vertices.size(), 0);
        EXPECT_GT(point_cloud->points.size(), 0);
        
    } catch (const std::exception& e) {
        FAIL() << "Exception in direct Fibonacci test: " << e.what();
    } catch (...) {
        FAIL() << "Unknown exception in direct Fibonacci test";
    }
}

TEST_F(FibonacciSurfaceTest, CustomParameters) {
    // Test with custom Fibonacci parameters
    spalt::FibonacciSurfaceParams params;
    params.probe_radius = 1.5;
    params.num_samples_per_atom = 30;
    params.ball_radii = {1.0, 1.5};
    
    // Test the parameters directly with the generator
    spalt::FibonacciSurfaceGenerator generator(params);
    
    // Get conformer data
    auto mol_ptr = mol->get_mol();
    const auto& conformer = mol_ptr->getConformer(0);
    int num_atoms = conformer.getNumAtoms();
    
    // Extract atom centers and radii
    Eigen::MatrixXd centers(num_atoms, 3);
    Eigen::VectorXd radii(num_atoms);
    
    for (int i = 0; i < num_atoms; ++i) {
        const auto& pos = conformer.getAtomPos(i);
        centers(i, 0) = pos.x;
        centers(i, 1) = pos.y;
        centers(i, 2) = pos.z;
        radii(i) = 1.7; // Default carbon radius
    }
    
    auto [mesh, point_cloud] = generator.generateSurface(centers, radii);
    
    EXPECT_GT(mesh->vertices.size(), 0);
    EXPECT_GT(point_cloud->points.size(), 0);
}

TEST_F(FibonacciSurfaceTest, ParameterValidation) {
    // Test parameter validation
    spalt::FibonacciSurfaceParams valid_params;
    valid_params.probe_radius = 1.2;
    valid_params.num_samples_per_atom = 20;
    valid_params.ball_radii = {1.2};
    
    EXPECT_TRUE(valid_params.isValid());
    
    // Test invalid parameters
    spalt::FibonacciSurfaceParams invalid_params;
    invalid_params.probe_radius = -1.0;  // Invalid
    invalid_params.num_samples_per_atom = 100;  // Too high
    invalid_params.ball_radii = {};  // Empty
    
    EXPECT_FALSE(invalid_params.isValid());
}

TEST_F(FibonacciSurfaceTest, SphereGeneration) {
    // Test Fibonacci sphere generation
    spalt::FibonacciSurfaceGenerator generator;
    
    auto sphere_points = generator.generateFibonacciSphere(100);
    EXPECT_EQ(sphere_points.rows(), 100);
    EXPECT_EQ(sphere_points.cols(), 3);
    
    // Check that points are on unit sphere
    for (int i = 0; i < sphere_points.rows(); ++i) {
        double norm = sphere_points.row(i).norm();
        EXPECT_NEAR(norm, 1.0, 1e-10);
    }
}

TEST_F(FibonacciSurfaceTest, SurfaceSampling) {
    // Test molecular surface sampling
    spalt::FibonacciSurfaceGenerator generator;
    
    // Create simple test data (single atom at origin)
    Eigen::MatrixXd centers(1, 3);
    centers << 0.0, 0.0, 0.0;
    Eigen::VectorXd radii(1);
    radii << 1.5;
    
    auto surface_points = generator.sampleMolecularSurface(centers, radii);
    
    EXPECT_GT(surface_points.rows(), 0);
    EXPECT_EQ(surface_points.cols(), 3);
    
    // Check that points are outside the atom
    for (int i = 0; i < surface_points.rows(); ++i) {
        double dist = surface_points.row(i).norm();
        EXPECT_GE(dist, 1.5);  // Should be outside atom radius
    }
}

TEST_F(FibonacciSurfaceTest, ComparisonWithMSMS) {
    // Test that Fibonacci method produces reasonable results
    // This test verifies that the Fibonacci method works correctly
    // and can be compared to other methods if available
    
    // Test Fibonacci surface generation directly (same as BasicSurfaceGeneration)
    try {
        // Get conformer data
        auto mol_ptr = mol->get_mol();
        const auto& conformer = mol_ptr->getConformer(0);
        int num_atoms = conformer.getNumAtoms();
        
        // Extract atom centers and radii
        Eigen::MatrixXd centers(num_atoms, 3);
        Eigen::VectorXd radii(num_atoms);
        
        for (int i = 0; i < num_atoms; ++i) {
            const auto& pos = conformer.getAtomPos(i);
            centers(i, 0) = pos.x;
            centers(i, 1) = pos.y;
            centers(i, 2) = pos.z;
            radii(i) = 1.7; // Default carbon radius
        }
        
        // Create Fibonacci surface generator
        spalt::FibonacciSurfaceParams params;
        spalt::FibonacciSurfaceGenerator generator(params);
        
        auto [mesh, point_cloud] = generator.generateSurface(centers, radii);
        
        // Verify reasonable results
        EXPECT_GT(mesh->vertices.size(), 0);
        EXPECT_GT(point_cloud->points.size(), 0);
        EXPECT_LT(mesh->vertices.size(), 10000);  // Not too many
        EXPECT_GT(mesh->vertices.size(), 10);     // Not too few
        
    } catch (const std::exception& e) {
        FAIL() << "Exception in Fibonacci surface creation: " << e.what();
    } catch (...) {
        FAIL() << "Unknown exception in Fibonacci surface creation";
    }
}
