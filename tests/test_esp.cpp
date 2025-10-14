#include <gtest/gtest.h>
#include <filesystem>
#include "props/esp.h"
#include "surface.h"

class ESPTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Clean up any leftover temp directories
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(".")) {
            if (entry.is_directory() && entry.path().filename().string().find("temp_esp_") == 0) {
                std::filesystem::remove_all(entry.path(), ec);
            }
        }
    }

    void TearDown() override {
        // Clean up any remaining temp directories
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(".")) {
            if (entry.is_directory() && entry.path().filename().string().find("temp_esp_") == 0) {
                std::filesystem::remove_all(entry.path(), ec);
            }
        }
    }
};

TEST_F(ESPTest, Constructor) {
    SurfaceESPProperty esp_prop;
    EXPECT_EQ(esp_prop.key(), "surface.esp");
}

TEST_F(ESPTest, LoadESPFromFile) {
    SurfaceESPProperty esp_prop;

    // Create a temporary ESP file for testing
    std::string test_file = "test_esp_data.txt";
    std::ofstream file(test_file);
    file << "0.1\n";
    file << "0.2\n";
    file << "-0.3\n";
    file.close();

    try {
        auto charges = esp_prop.load_esp(test_file);
        EXPECT_EQ(charges.size(), 3);
        EXPECT_DOUBLE_EQ(charges[0], 0.1);
        EXPECT_DOUBLE_EQ(charges[1], 0.2);
        EXPECT_DOUBLE_EQ(charges[2], -0.3);

        // Clean up
        std::filesystem::remove(test_file);
    } catch (const std::exception& e) {
        std::filesystem::remove(test_file);
        FAIL() << "Exception: " << e.what();
    }
}

TEST_F(ESPTest, LoadESPFromNonexistentFile) {
    SurfaceESPProperty esp_prop;

    EXPECT_THROW(esp_prop.load_esp("nonexistent_file.txt"), std::runtime_error);
}

TEST_F(ESPTest, GetESPChargesFromCache) {
    std::unordered_map<std::string, std::any> cache;
    std::vector<double> test_charges = {0.1, 0.2, -0.3};
    cache["esp_charges"] = test_charges;

    auto retrieved_charges = SurfaceESPProperty::get_esp_charges(cache);
    EXPECT_EQ(retrieved_charges.size(), 3);
    EXPECT_DOUBLE_EQ(retrieved_charges[0], 0.1);
    EXPECT_DOUBLE_EQ(retrieved_charges[1], 0.2);
    EXPECT_DOUBLE_EQ(retrieved_charges[2], -0.3);
}

TEST_F(ESPTest, GetESPChargesFromEmptyCache) {
    std::unordered_map<std::string, std::any> cache;

    EXPECT_THROW(SurfaceESPProperty::get_esp_charges(cache), std::runtime_error);
}

TEST_F(ESPTest, ConstructorBasedConfiguration) {
    // Test the new extensible API with constructor-based configuration
    Molecule molecule("CCO");  // Ethanol

    // Create ESP property with RDKit method using the new API
    SurfaceESPProperty esp_prop(ChargeMethod::RDKIT);
    EXPECT_EQ(esp_prop.getChargeMethod(), ChargeMethod::RDKIT);

    // Register the property and compute using new API
    molecule.registerProperty(
        "esp_rdkit", []() { return std::make_unique<SurfaceESPProperty>(ChargeMethod::RDKIT); });

    // Generate conformers first
    molecule.calculate_conformers(1);
    EXPECT_NO_THROW(molecule.compute("esp_rdkit", 0));

    // Test that we can get the results
    auto charges = molecule.getProperty<std::vector<double>>("esp_rdkit", "esp_charges");
    auto potentials =
        molecule.getProperty<std::vector<double>>("esp_rdkit", "electrostatic_potential");

    EXPECT_GT(charges.size(), 0);
    EXPECT_GT(potentials.size(), 0);

    std::cout << "Constructor-based ESP computation successful: " << charges.size() << " charges, "
              << potentials.size() << " potentials" << std::endl;
}

TEST_F(ESPTest, DefaultConstructorUsesXTB) {
    SurfaceESPProperty esp_prop;
    EXPECT_EQ(esp_prop.getChargeMethod(), ChargeMethod::XTB);
}

TEST_F(ESPTest, ComputeDefaultMethod) {
    // Test the new extensible API
    Molecule molecule("CCO");  // Ethanol

    // Register ESP property with XTB method (might fail if not installed)
    molecule.registerProperty(
        "esp_xtb", []() { return std::make_unique<SurfaceESPProperty>(ChargeMethod::XTB); });

    // This will try to run XTB, which might fail if not installed
    // We expect it to either succeed or throw a runtime_error
    try {
        molecule.compute("esp_xtb", 0);
        // If successful, check that we can get the results
        auto charges = molecule.getProperty<std::vector<double>>("esp_xtb", "esp_charges");
        auto potentials =
            molecule.getProperty<std::vector<double>>("esp_xtb", "electrostatic_potential");
        EXPECT_GT(charges.size(), 0);
        EXPECT_GT(potentials.size(), 0);
        std::cout << "XTB ESP computation successful: " << charges.size() << " charges, "
                  << potentials.size() << " potentials" << std::endl;
    } catch (const std::runtime_error& e) {
        // Expected if xtb is not installed
        std::cout << "XTB computation failed (expected if xtb not installed): " << e.what()
                  << std::endl;
    }
}

TEST_F(ESPTest, RDKitChargeComputation) {
    // Test RDKit charge computation using the new extensible API
    Molecule molecule("CCO");  // Ethanol

    // Generate conformers first
    molecule.calculate_conformers(1);
    // Use the default ESP property (which uses RDKit)
    EXPECT_NO_THROW(molecule.compute("esp", 0));

    // Test that we can get the results
    auto charges = molecule.getProperty<std::vector<double>>("esp", "esp_charges");
    auto potentials = molecule.getProperty<std::vector<double>>("esp", "electrostatic_potential");

    EXPECT_GT(charges.size(), 0);
    EXPECT_GT(potentials.size(), 0);

    std::cout << "RDKit ESP computation successful: " << charges.size() << " charges, "
              << potentials.size() << " potentials" << std::endl;
}

TEST_F(ESPTest, ConvenienceFunctions) {
    // Test the convenience functions for property registration
    Molecule molecule("CCO");  // Ethanol

    // Test getRegisteredProperties before registration
    auto properties_before = molecule.getRegisteredProperties();
    EXPECT_GT(properties_before.size(), 0);  // Should have default properties

    // Test registerAllProperties
    molecule.registerAllProperties();
    auto properties_after = molecule.getRegisteredProperties();
    EXPECT_GE(properties_after.size(), properties_before.size());

    // Generate conformers first
    molecule.calculate_conformers(1);
    // Should be able to compute with registered properties
    EXPECT_NO_THROW(molecule.compute("esp", 0));

    // Test that we can get the results
    auto charges = molecule.getProperty<std::vector<double>>("esp", "esp_charges");
    EXPECT_GT(charges.size(), 0);

    std::cout << "Convenience functions test successful. Registered properties: ";
    for (const auto& prop : properties_after) {
        std::cout << prop << " ";
    }
    std::cout << std::endl;
}

// Note: We don't test XTB computation here because it requires xtb to be installed
// and would create temporary directories. That should be tested in integration tests.
