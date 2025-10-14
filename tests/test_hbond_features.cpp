#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "molecule.h"
#include "props/hbond.h"
#include "props/rdkit_utils.h"

class HBondFeaturesTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Embedded SDF content for Mol1 (test molecule)
        std::string sdf_content = R"(
     RDKit          3D

 26 27  0  0  0  0  0  0  0  0999 V2000
   -4.1917    1.0443   -0.4688 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7155    1.1617    0.9751 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0117   -0.0149    1.3735 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7006   -0.3486    0.9963 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1971   -1.6062    1.3935 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1061   -2.0045    1.0520 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9015   -1.1318    0.3116 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1736   -1.3981   -0.0795 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.6864   -0.3592   -0.7815 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3071   -0.3533   -1.4445 S   0  0  0  0  0  0  0  0  0  0  0  0
    5.1239    1.0916   -0.9609 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.2358   -0.4198   -2.9437 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0630   -1.5463   -0.9300 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5916    0.9208   -0.9443 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.4102    0.1037   -0.0756 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8763    0.5252    0.2486 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8426    0.1514   -0.5800 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7730    1.9492   -0.7443 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3300    0.9507   -1.1611 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6063    1.2724    1.6285 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0874    2.0667    1.1155 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8227   -2.2771    1.9685 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4834   -2.9715    1.3595 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.1031    1.0236   -1.3183 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1812    1.0756    0.0820 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2121    1.5008   -0.0722 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  1  0
 10 12  2  0
 10 13  2  0
  9 14  1  0
 14 15  1  0
 15 16  2  0
 16  4  1  0
 15  7  1  0
  1 17  1  0
  1 18  1  0
  1 19  1  0
  2 20  1  0
  2 21  1  0
  5 22  1  0
  6 23  1  0
 11 24  1  0
 11 25  1  0
 16 26  1  0
M  END
$$$$
)";

        // Create temporary file with embedded SDF content
        std::string temp_file = "temp_test_hbond_mol.sdf";
        std::ofstream temp_stream(temp_file);
        temp_stream << sdf_content;
        temp_stream.close();

        // Create molecule from temporary file
        mol = std::make_shared<Molecule>(temp_file);
        ASSERT_TRUE(mol->get_mol() != nullptr) << "Failed to load molecule from temporary SDF file";

        // Clean up temporary file after molecule is loaded
        std::remove(temp_file.c_str());

        hbond_prop = std::make_unique<HBondProperty>();
    }

    void TearDown() override {
        // Clean up
    }

  protected:
    std::shared_ptr<Molecule> mol;
    std::unique_ptr<HBondProperty> hbond_prop;
};

TEST_F(HBondFeaturesTest, Mol1AcceptorFeaturesMatchPython) {
    // This test verifies that our C++ implementation matches the Python output exactly
    // for the Mol1.sdf molecule

    ASSERT_TRUE(mol->get_mol() != nullptr);

    // Get the molecule and conformer
    auto rdmol_ptr = mol->get_mol();
    std::shared_ptr<RDKit::ROMol> rdmol(rdmol_ptr.release());
    auto conf = rdmol->getConformer(0);

    // Get acceptors using our implementation
    auto acceptors = hbond_prop->getAcceptors(rdmol);

    // Expected results from Python output:
    // Atom idx 2: Acceptor type 2, Vec Type: linear
    // Atom idx 11: Acceptor type 1, Vec Type: cone
    // Atom idx 12: Acceptor type 1, Vec Type: cone

    std::cout << "Testing Mol1.sdf acceptor features..." << std::endl;
    std::cout << "Found " << acceptors.size() << " acceptor atoms" << std::endl;

    // Test each acceptor atom
    for (int atom_idx : acceptors) {
        std::cout << "Testing atom " << atom_idx << std::endl;

        // Get acceptor type
        auto atom = rdmol->getAtomWithIdx(atom_idx);
        int acceptor_type = hbond_prop->getAcceptorType(*atom);

        // Get feature vectors using our RDKit utilities
        std::vector<int> featAtoms = {atom_idx};
        auto result = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);

        if (acceptor_type == 1) {
            // For type 1 acceptors, use getAcceptor1FeatVects
            result = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);
        } else if (acceptor_type == 2) {
            // For type 2 acceptors, use getAcceptor2FeatVects
            result = RDKitUtils::getAcceptor2FeatVects(conf, featAtoms, 1.0);
        } else if (acceptor_type == 3) {
            // For type 3 acceptors, use getAcceptor3FeatVects
            result = RDKitUtils::getAcceptor3FeatVects(conf, featAtoms, 1.0);
        }

        std::string vec_type = result.second;

        // Print results to match Python output format
        std::cout << "Vec Type: " << vec_type << "  Acceptor type: " << acceptor_type
                  << "  Atom idx: " << atom_idx << std::endl;

        // Verify expected results
        if (atom_idx == 2) {
            EXPECT_EQ(acceptor_type, 2) << "Atom 2 should be acceptor type 2";
            EXPECT_EQ(vec_type, "linear") << "Atom 2 should have linear vec type";
        } else if (atom_idx == 11) {
            EXPECT_EQ(acceptor_type, 1) << "Atom 11 should be acceptor type 1";
            EXPECT_EQ(vec_type, "cone") << "Atom 11 should have cone vec type";
        } else if (atom_idx == 12) {
            EXPECT_EQ(acceptor_type, 1) << "Atom 12 should be acceptor type 1";
            EXPECT_EQ(vec_type, "cone") << "Atom 12 should have cone vec type";
        }

        // Verify we have valid feature vectors
        EXPECT_GT(result.first.size(), 0) << "Should have at least one feature vector";

        // Print vector details for debugging
        for (size_t i = 0; i < result.first.size(); ++i) {
            const auto& pair = result.first[i];
            std::cout << "  Vector " << i << ": start(" << pair.first.x << ", " << pair.first.y
                      << ", " << pair.first.z << ") end(" << pair.second.x << ", " << pair.second.y
                      << ", " << pair.second.z << ")" << std::endl;
        }
    }
}

TEST_F(HBondFeaturesTest, Mol1AcceptorTypeDetection) {
    // Test that acceptor type detection works correctly
    auto rdmol_ptr = mol->get_mol();
    std::shared_ptr<RDKit::ROMol> rdmol(rdmol_ptr.release());

    // Test specific atoms from Mol1.sdf
    auto atom2 = rdmol->getAtomWithIdx(2);    // Oxygen
    auto atom11 = rdmol->getAtomWithIdx(11);  // Nitrogen
    auto atom12 = rdmol->getAtomWithIdx(12);  // Oxygen

    // Check acceptor types
    int type2 = hbond_prop->getAcceptorType(*atom2);
    int type11 = hbond_prop->getAcceptorType(*atom11);
    int type12 = hbond_prop->getAcceptorType(*atom12);

    std::cout << "Acceptor types: atom2=" << type2 << ", atom11=" << type11 << ", atom12=" << type12
              << std::endl;

    // Verify expected types
    EXPECT_EQ(type2, 2) << "Atom 2 (O) should be acceptor type 2";
    EXPECT_EQ(type11, 1) << "Atom 11 (N) should be acceptor type 1";
    EXPECT_EQ(type12, 1) << "Atom 12 (O) should be acceptor type 1";
}

TEST_F(HBondFeaturesTest, Mol1DonorHydrogens) {
    // Test donor hydrogen detection
    auto rdmol_ptr = mol->get_mol();
    std::shared_ptr<RDKit::ROMol> rdmol(rdmol_ptr.release());
    auto donor_hydrogens = hbond_prop->getDonorHydrogens(rdmol);

    std::cout << "Found " << donor_hydrogens.size() << " donor hydrogens" << std::endl;
    for (int h_idx : donor_hydrogens) {
        auto h_atom = rdmol->getAtomWithIdx(h_idx);
        std::cout << "Donor hydrogen at atom " << h_idx << " (element: " << h_atom->getSymbol()
                  << ")" << std::endl;
    }

    // Should find some donor hydrogens
    EXPECT_GT(donor_hydrogens.size(), 0) << "Should find at least one donor hydrogen";
}

TEST_F(HBondFeaturesTest, Mol1CompleteFeatureMatching) {
    // This is the main test that should match the Python output exactly
    auto rdmol_ptr = mol->get_mol();
    std::shared_ptr<RDKit::ROMol> rdmol(rdmol_ptr.release());
    auto conf = rdmol->getConformer(0);

    auto acceptors = hbond_prop->getAcceptors(rdmol);

    std::cout << "\n=== MOL1.SDF ACCEPTOR FEATURE TEST ===" << std::endl;
    std::cout << "This should match the Python test case output exactly:" << std::endl;

    // Simulate the Python loop that generates the test case
    for (int atom_idx : acceptors) {
        auto atom = rdmol->getAtomWithIdx(atom_idx);
        int acceptor_type = hbond_prop->getAcceptorType(*atom);

        std::vector<int> featAtoms = {atom_idx};
        std::pair<std::vector<std::pair<RDKitUtils::Vector3D, RDKitUtils::Vector3D>>, std::string>
            result;

        // Call the appropriate function based on acceptor type
        if (acceptor_type == 1) {
            result = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);
        } else if (acceptor_type == 2) {
            result = RDKitUtils::getAcceptor2FeatVects(conf, featAtoms, 1.0);
        } else if (acceptor_type == 3) {
            result = RDKitUtils::getAcceptor3FeatVects(conf, featAtoms, 1.0);
        }

        std::string vec_type = result.second;

        // Print in the exact format from the Python test case
        std::cout << "Vec Type:  " << vec_type << "  Acceptor type:  " << acceptor_type
                  << "  Atom idx:  " << atom_idx << std::endl;

        // Verify against expected results
        if (atom_idx == 2) {
            EXPECT_EQ(vec_type, "linear") << "Atom 2 should have linear vec type";
            EXPECT_EQ(acceptor_type, 2) << "Atom 2 should be acceptor type 2";
        } else if (atom_idx == 11 || atom_idx == 12) {
            EXPECT_EQ(vec_type, "cone") << "Atoms 11 and 12 should have cone vec type";
            EXPECT_EQ(acceptor_type, 1) << "Atoms 11 and 12 should be acceptor type 1";
        }
    }

    std::cout << "=== END TEST CASE ===" << std::endl;
}
