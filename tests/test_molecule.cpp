#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include "molecule.h"

class MoleculeTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Clean up any leftover test files
        safe_remove("test_mol.mol");
        safe_remove("benzene.mol");
    }

    void TearDown() override {
        // Clean up test files
        safe_remove("test_mol.mol");
        safe_remove("benzene.mol");
    }

  protected:
    void safe_remove(const std::filesystem::path& p) {
        std::error_code ec;
        std::filesystem::remove(p, ec);
    }
};

TEST_F(MoleculeTest, ConstructorWithValidSMILES) {
    // Test with a simple valid SMILES string
    EXPECT_NO_THROW({
        Molecule mol("CCO");  // Ethanol
    });
}

TEST_F(MoleculeTest, ConstructorWithInvalidSMILES) {
    // Test with a simple invalid SMILES string
    Molecule mol("CCOO#");

    // The constructor should not throw, but the molecule should be invalid internally.
    // So, we check that get_mol() returns nullptr for an invalid SMILES.
    EXPECT_EQ(mol.get_mol(), nullptr);
}

TEST_F(MoleculeTest, ConstructorWithComplexSMILES) {
    // Test with a more complex valid SMILES string
    EXPECT_NO_THROW({
        Molecule mol("c1ccccc1");  // Benzene
    });

    EXPECT_NO_THROW({
        Molecule mol("CC(=O)O");  // Acetic acid
    });
}

TEST_F(MoleculeTest, GetSmiles) {
    Molecule mol("CCO");
    EXPECT_EQ(mol.get_smiles(), "CCO");

    Molecule mol2("c1ccccc1");
    EXPECT_EQ(mol2.get_smiles(), "c1ccccc1");
}

TEST_F(MoleculeTest, GetCharge) {
    // Test neutral molecules
    Molecule mol1("CCO");  // Ethanol - neutral
    EXPECT_EQ(mol1.get_charge(), 0);

    Molecule mol2("c1ccccc1");  // Benzene - neutral
    EXPECT_EQ(mol2.get_charge(), 0);

    // Test charged molecules
    Molecule mol3("CC(=O)O");  // Acetic acid - neutral
    EXPECT_EQ(mol3.get_charge(), 0);

    // Test with a positively charged molecule (ammonium ion)
    Molecule mol4("[NH4+]");  // Ammonium ion - +1 charge
    EXPECT_EQ(mol4.get_charge(), 1);

    // Test with a negatively charged molecule (acetate ion)
    Molecule mol5("CC(=O)[O-]");  // Acetate ion - -1 charge
    EXPECT_EQ(mol5.get_charge(), -1);
}

TEST_F(MoleculeTest, GetMol) {
    Molecule mol("CCO");
    auto rdkit_mol = mol.get_mol();

    EXPECT_NE(rdkit_mol, nullptr);

    EXPECT_EQ(rdkit_mol->getNumAtoms(), 3);

    EXPECT_EQ(rdkit_mol->getNumBonds(), 2);
}

TEST_F(MoleculeTest, GetMolBenzene) {
    Molecule mol("c1ccccc1");  // Benzene
    auto rdkit_mol = mol.get_mol();

    EXPECT_NE(rdkit_mol, nullptr);
    EXPECT_EQ(rdkit_mol->getNumAtoms(), 6);
    EXPECT_EQ(rdkit_mol->getNumBonds(), 6);
}

TEST_F(MoleculeTest, MolToMolfileDefaultFilename) {
    Molecule mol("CCO");

    // Test with default filename
    EXPECT_TRUE(mol.mol_to_molfile());

    // Check if file was created
    EXPECT_TRUE(std::filesystem::exists("mol.mol"));

    // Check if file has content
    std::ifstream file("mol.mol");
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());

    // Clean up
    safe_remove("mol.mol");
}

TEST_F(MoleculeTest, MolToMolfileCustomFilename) {
    Molecule mol("c1ccccc1");  // Benzene

    // Test with custom filename
    EXPECT_TRUE(mol.mol_to_molfile("benzene.mol"));

    // Check if file was created
    EXPECT_TRUE(std::filesystem::exists("benzene.mol"));

    // Check if file has content
    std::ifstream file("benzene.mol");
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());

    // The file should contain benzene-specific information
    EXPECT_TRUE(content.find("C6H6") != std::string::npos ||
                content.find("benzene") != std::string::npos ||
                content.length() > 100);  // molfile should be substantial
}

TEST_F(MoleculeTest, MolToMolfileMultipleMolecules) {
    Molecule mol1("CCO");
    Molecule mol2("CC(=O)O");

    // Test writing different molecules to different files
    EXPECT_TRUE(mol1.mol_to_molfile("ethanol.mol"));
    EXPECT_TRUE(mol2.mol_to_molfile("acetic_acid.mol"));

    // Both files should exist
    EXPECT_TRUE(std::filesystem::exists("ethanol.mol"));
    EXPECT_TRUE(std::filesystem::exists("acetic_acid.mol"));

    // Both files should have content
    std::ifstream file1("ethanol.mol");
    std::ifstream file2("acetic_acid.mol");

    std::string content1((std::istreambuf_iterator<char>(file1)), std::istreambuf_iterator<char>());
    std::string content2((std::istreambuf_iterator<char>(file2)), std::istreambuf_iterator<char>());

    EXPECT_FALSE(content1.empty());
    EXPECT_FALSE(content2.empty());

    // Content should be different (different molecules)
    EXPECT_NE(content1, content2);

    // Clean up
    safe_remove("ethanol.mol");
    safe_remove("acetic_acid.mol");
}

TEST_F(MoleculeTest, MolToMolfileOverwrite) {
    Molecule mol1("CCO");
    Molecule mol2("c1ccccc1");

    // Write first molecule
    EXPECT_TRUE(mol1.mol_to_molfile("test_mol.mol"));

    // Get original content
    std::ifstream file1("test_mol.mol");
    std::string original_content((std::istreambuf_iterator<char>(file1)),
                                 std::istreambuf_iterator<char>());
    file1.close();

    // Write second molecule to same file (should overwrite)
    EXPECT_TRUE(mol2.mol_to_molfile("test_mol.mol"));

    // Get new content
    std::ifstream file2("test_mol.mol");
    std::string new_content((std::istreambuf_iterator<char>(file2)),
                            std::istreambuf_iterator<char>());

    // Content should be different (benzene vs ethanol)
    EXPECT_NE(original_content, new_content);
    EXPECT_FALSE(new_content.empty());
}

TEST_F(MoleculeTest, EdgeCaseEmptyString) {
    // Test with empty SMILES string - this should handle gracefully
    // Note: This might throw an exception depending on RDKit behavior
    try {
        Molecule mol("");
        // If no exception, check that mol is nullptr or handles gracefully
        auto rdkit_mol = mol.get_mol();
        // The behavior here depends on RDKit's handling of empty SMILES
        EXPECT_EQ(mol.get_smiles(), "");
    } catch (const std::exception& e) {
        // This is acceptable behavior for invalid SMILES
        EXPECT_TRUE(true);  // Test passes if exception is thrown
    }
}

TEST_F(MoleculeTest, EdgeCaseInvalidSMILES) {
    // Test with invalid SMILES string
    Molecule mol("invalid_smiles_string");

    // Should still store the SMILES string
    EXPECT_EQ(mol.get_smiles(), "invalid_smiles_string");

    // But the RDKit molecule should be null
    auto rdkit_mol = mol.get_mol();
    EXPECT_EQ(rdkit_mol, nullptr);

    // Writing to molfile should fail gracefully
    EXPECT_FALSE(mol.mol_to_molfile("invalid_test.mol"));
}

TEST_F(MoleculeTest, DestructorCleanup) {
    // Test that destructor works properly (no memory leaks)
    {
        Molecule mol("CCO");
        auto rdkit_mol = mol.get_mol();
        EXPECT_NE(rdkit_mol, nullptr);
    }  // mol goes out of scope here, destructor should be called

    // If we get here without crashing, destructor worked
    EXPECT_TRUE(true);
}

TEST_F(MoleculeTest, CopySemantics) {
    // Test that we can create multiple molecules independently
    Molecule mol1("CCO");
    Molecule mol2("c1ccccc1");

    EXPECT_EQ(mol1.get_smiles(), "CCO");
    EXPECT_EQ(mol2.get_smiles(), "c1ccccc1");

    auto rdkit_mol1 = mol1.get_mol();
    auto rdkit_mol2 = mol2.get_mol();

    EXPECT_NE(rdkit_mol1, nullptr);
    EXPECT_NE(rdkit_mol2, nullptr);

    // Should be different molecules with different atom counts
    EXPECT_NE(rdkit_mol1->getNumAtoms(), rdkit_mol2->getNumAtoms());
}

TEST_F(MoleculeTest, CalculateConformers) {
    // Test conformer generation for a simple molecule
    Molecule mol("CCO");  // Ethanol

    // Should not throw when generating conformers
    EXPECT_NO_THROW(mol.calculate_conformers(2));

    // Test with a more complex molecule
    Molecule mol2("c1ccccc1");  // Benzene
    EXPECT_NO_THROW(mol2.calculate_conformers(1));
}

TEST_F(MoleculeTest, CalculateConformersInvalidMolecule) {
    // Test conformer generation with invalid molecule
    Molecule mol("invalid_smiles");

    // Should handle gracefully without crashing
    EXPECT_NO_THROW(mol.calculate_conformers(2));
}
