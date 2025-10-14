#include <gtest/gtest.h>
#include "molecule.h"
#include "surface.h"
#include "conformer_generator.h"

class CLIIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create molecule like CLI does: SMILES + generate conformers
        test_molecule = std::make_unique<Molecule>("C1CCCCC1");
        // Generate conformers like CLI does
        if (test_molecule->get_mol() != nullptr) {
            test_molecule->calculate_conformers(5);  // Same as CLI default
        }
    }

    std::unique_ptr<Molecule> test_molecule;
};

TEST_F(CLIIntegrationTest, FibonacciSurfaceCreationLikeCLI) {
    // Test the exact sequence that CLI uses
    ASSERT_NE(test_molecule, nullptr);
    EXPECT_GT(test_molecule->getNumConformers(), 0);
    
    // This mimics exactly what CLI does: createSurface with fibonacci
    EXPECT_NO_THROW({
        int surface_id = test_molecule->createSurface(
            0,          // conformer_id
            1000,       // num_vertices 
            1.2,        // radius
            3.0,        // density
            3.0,        // hdensity
            "tses",     // type
            "fibonacci", // program
            "full"      // sample_method
        );
        EXPECT_GE(surface_id, 0);
    });
}

TEST_F(CLIIntegrationTest, DirectSurfaceConstructorWithFibonacci) {
    // Test direct Surface constructor with fibonacci (what createSurface calls internally)
    ASSERT_NE(test_molecule, nullptr);
    EXPECT_GT(test_molecule->getNumConformers(), 0);
    
    EXPECT_NO_THROW({
        Surface surface(*test_molecule, 0, 1000, 1.2, 3.0, 3.0, "tses", "fibonacci", "full");
        // If we get here, the surface was created successfully
        SUCCEED();
    });
}

TEST_F(CLIIntegrationTest, CompareWithWorkingMSMS) {
    // Test that MSMS works (to establish baseline)
    ASSERT_NE(test_molecule, nullptr);
    EXPECT_GT(test_molecule->getNumConformers(), 0);
    
    // MSMS should work
    EXPECT_NO_THROW({
        int surface_id = test_molecule->createSurface(
            0, 1000, 1.2, 3.0, 3.0, "tses", "msms", "full"
        );
        EXPECT_GE(surface_id, 0);
    });
}

