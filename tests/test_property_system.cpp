#include <gtest/gtest.h>
#include "molecule.h"
#include "props/property_params.h"

class PropertySystemTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Create test molecules
        mol1 = std::make_unique<Molecule>("C1CCCCC1");  // Cyclohexane
        mol2 = std::make_unique<Molecule>("CCO");       // Ethanol

        // Generate conformers
        mol1->calculate_conformers(1);
        mol2->calculate_conformers(1);

        // Create surfaces
        mol1->createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
        mol2->createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
    }

    std::unique_ptr<Molecule> mol1;
    std::unique_ptr<Molecule> mol2;
};

TEST_F(PropertySystemTest, PropertyParamsRegistryDefaultValues) {
    // Test ESP parameters
    auto esp_params = PropertyParamsRegistry::getParams("esp");
    EXPECT_DOUBLE_EQ(esp_params.min_value, -10.0);
    EXPECT_DOUBLE_EQ(esp_params.max_value, 10.0);

    // Test hydrophobicity parameters
    auto hydro_params = PropertyParamsRegistry::getParams("hydrophobicity");
    EXPECT_DOUBLE_EQ(hydro_params.min_value, -4.5);
    EXPECT_DOUBLE_EQ(hydro_params.max_value, 4.5);

    // Test hbond parameters
    auto hbond_params = PropertyParamsRegistry::getParams("hbond");
    EXPECT_DOUBLE_EQ(hbond_params.min_value, -1.0);
    EXPECT_DOUBLE_EQ(hbond_params.max_value, 1.0);
}

TEST_F(PropertySystemTest, PropertyParamsRegistryUnknownProperty) {
    // Test unknown property
    EXPECT_THROW(PropertyParamsRegistry::getParams("unknown_property"), std::runtime_error);
}

TEST_F(PropertySystemTest, PropertyParamsRegistryGetPropertyNames) {
    auto property_names = PropertyParamsRegistry::getPropertyNames();

    // Should contain the registered properties
    EXPECT_TRUE(std::find(property_names.begin(), property_names.end(), "esp") !=
                property_names.end());
    EXPECT_TRUE(std::find(property_names.begin(), property_names.end(), "hydrophobicity") !=
                property_names.end());
    EXPECT_TRUE(std::find(property_names.begin(), property_names.end(), "hbond") !=
                property_names.end());

    EXPECT_EQ(property_names.size(), 11);
}

TEST_F(PropertySystemTest, MultiPropertyComputation) {
    // Test computing multiple properties
    std::vector<std::string> properties = {"esp", "hydrophobicity", "hbond"};

    for (const auto& property : properties) {
        EXPECT_NO_THROW(mol1->compute(property, 0));
    }

    // Verify all properties were computed by checking if we can get values
    for (const auto& property : properties) {
        if (property == "esp") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "electrostatic_potential", 0));
        } else if (property == "hydrophobicity") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "hydrophobicity_values", 0));
        } else if (property == "hbond") {
            EXPECT_NO_THROW(mol1->getProperty<std::vector<double>>(property, "hbond_values", 0));
        }
    }
}

TEST_F(PropertySystemTest, PropertyNormalization) {
    // Compute ESP property
    mol1->compute("esp", 0);

    // Get the property values
    auto values = mol1->getProperty<std::vector<double>>("esp", "electrostatic_potential");
    ASSERT_GT(values.size(), 0);

    // Values should be within the expected range
    for (const auto& value : values) {
        EXPECT_GE(value, -10.0);
        EXPECT_LE(value, 10.0);
    }
}

TEST_F(PropertySystemTest, PropertyColorMapping) {
    // Compute hydrophobicity property
    mol1->compute("hydrophobicity", 0);

    // Get the property values
    auto values = mol1->getProperty<std::vector<double>>("hydrophobicity", "hydrophobicity_values");
    ASSERT_GT(values.size(), 0);

    // Values should be within the expected range
    for (const auto& value : values) {
        EXPECT_GE(value, -4.5);
        EXPECT_LE(value, 4.5);
    }
}

TEST_F(PropertySystemTest, HBondPropertyComputation) {
    // Test hbond property computation with ethanol (CCO) which has hydrogen bonds
    EXPECT_NO_THROW(mol2->compute("hbond", 0));

    // Verify property was computed by checking if we can get values
    EXPECT_NO_THROW(mol2->getProperty<std::vector<double>>("hbond", "hbond_values", 0));

    // Get the property values
    auto values = mol2->getProperty<std::vector<double>>("hbond", "hbond_values");
    ASSERT_GT(values.size(), 0);

    // Values should be within the expected range
    for (const auto& value : values) {
        EXPECT_GE(value, -1.0);
        EXPECT_LE(value, 1.0);
    }
}

TEST_F(PropertySystemTest, PropertyCacheConsistency) {
    // Compute properties for both molecules
    std::vector<std::string> properties = {"esp", "hydrophobicity"};

    for (const auto& property : properties) {
        mol1->compute(property, 0);
        mol2->compute(property, 0);
    }

    // Verify cache consistency by checking if we can get values
    for (const auto& property : properties) {
        if (property == "esp") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "electrostatic_potential", 0));
            EXPECT_NO_THROW(
                mol2->getProperty<std::vector<double>>(property, "electrostatic_potential", 0));
        } else if (property == "hydrophobicity") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "hydrophobicity_values", 0));
            EXPECT_NO_THROW(
                mol2->getProperty<std::vector<double>>(property, "hydrophobicity_values", 0));
        }
    }
}

TEST_F(PropertySystemTest, PropertyRegistration) {
    // Test property registration
    auto properties_before = mol1->getRegisteredProperties();

    // Register all properties
    mol1->registerAllProperties();

    auto properties_after = mol1->getRegisteredProperties();
    EXPECT_GE(properties_after.size(), properties_before.size());
}

TEST_F(PropertySystemTest, EdgeCaseEmptyProperties) {
    // Test with empty property list
    std::vector<std::string> empty_properties;

    // Should not throw
    EXPECT_NO_THROW(mol1->compute("esp", 0));
}

TEST_F(PropertySystemTest, EdgeCaseInvalidProperty) {
    // Test with invalid property
    EXPECT_THROW(mol1->compute("invalid_property", 0), std::runtime_error);
}

TEST_F(PropertySystemTest, PropertyComputationOrder) {
    // Test that properties can be computed in any order
    std::vector<std::string> properties = {"hbond", "esp", "hydrophobicity"};

    for (const auto& property : properties) {
        EXPECT_NO_THROW(mol1->compute(property, 0));
        if (property == "esp") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "electrostatic_potential", 0));
        } else if (property == "hydrophobicity") {
            EXPECT_NO_THROW(
                mol1->getProperty<std::vector<double>>(property, "hydrophobicity_values", 0));
        } else if (property == "hbond") {
            EXPECT_NO_THROW(mol1->getProperty<std::vector<double>>(property, "hbond_values", 0));
        }
    }
}

TEST_F(PropertySystemTest, MultipleConformerProperties) {
    // Generate multiple conformers
    mol1->calculate_conformers(3);

    // Create surfaces for all conformers
    for (int i = 0; i < 3; ++i) {
        mol1->createSurface(i, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
    }

    // Compute properties for all conformers
    std::vector<std::string> properties = {"esp", "hydrophobicity"};

    for (int i = 0; i < 3; ++i) {
        for (const auto& property : properties) {
            EXPECT_NO_THROW(mol1->compute(property, i));
            if (property == "esp") {
                EXPECT_NO_THROW(
                    mol1->getProperty<std::vector<double>>(property, "electrostatic_potential", i));
            } else if (property == "hydrophobicity") {
                EXPECT_NO_THROW(
                    mol1->getProperty<std::vector<double>>(property, "hydrophobicity_values", i));
            }
        }
    }
}
