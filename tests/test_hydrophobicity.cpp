#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <GraphMol/FileParsers/MolSupplier.h>
#include "molecule.h"
#include "props/hydrophobicity.h"

TEST(HydrophobicityProperty, TestMol1SDF) {
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
    std::string temp_file = "temp_test_mol.sdf";
    std::ofstream temp_stream(temp_file);
    temp_stream << sdf_content;
    temp_stream.close();

    // Create molecule from temporary file
    Molecule mol(temp_file);
    ASSERT_TRUE(mol.get_mol() != nullptr) << "Failed to load molecule from temporary SDF file";

    // Clean up temporary file after molecule is loaded
    std::remove(temp_file.c_str());

    // Create hydrophobicity property
    HydrophobicityProperty hydrophobicity;

    // Test the molAssignKd method directly
    auto atom_kd = hydrophobicity.molAssignKd(mol.get_mol());

    // Print results for verification
    std::cout << "Hydrophobicity values for embedded test molecule:" << std::endl;
    for (const auto& pair : atom_kd) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
    }

    // Basic validation - should have values for all atoms
    EXPECT_GT(atom_kd.size(), 0) << "No hydrophobicity values computed";

    // Check that we have values for expected atoms (based on the SDF structure)
    // The molecule has 26 atoms (C, O, N, S, H)
    EXPECT_EQ(atom_kd.size(), 26) << "Expected 26 atoms, got " << atom_kd.size();

    // Check for specific atoms that should exist
    EXPECT_TRUE(atom_kd.find("C0") != atom_kd.end()) << "Atom C0 not found";
    EXPECT_TRUE(atom_kd.find("C1") != atom_kd.end()) << "Atom C1 not found";
    EXPECT_TRUE(atom_kd.find("O2") != atom_kd.end()) << "Atom O2 not found";
    EXPECT_TRUE(atom_kd.find("N7") != atom_kd.end()) << "Atom N7 not found";
    EXPECT_TRUE(atom_kd.find("S9") != atom_kd.end()) << "Atom S9 not found";

    // Values should be within reasonable range
    for (const auto& pair : atom_kd) {
        EXPECT_GE(pair.second, -4.5) << "Value too low for " << pair.first;
        EXPECT_LE(pair.second, 4.5) << "Value too high for " << pair.first;
    }
}

TEST(HydrophobicityProperty, TestKdFromLogP) {
    HydrophobicityProperty hydrophobicity;

    // Test the kdFromLogP function with known values
    double kd1 = hydrophobicity.kdFromLogP(0.0);   // Neutral LogP
    double kd2 = hydrophobicity.kdFromLogP(1.0);   // Positive LogP
    double kd3 = hydrophobicity.kdFromLogP(-1.0);  // Negative LogP

    // Values should be different
    EXPECT_NE(kd1, kd2) << "Different LogP values should give different kd values";
    EXPECT_NE(kd1, kd3) << "Different LogP values should give different kd values";
    EXPECT_NE(kd2, kd3) << "Different LogP values should give different kd values";

    // Values should be within expected range
    EXPECT_GE(kd1, -4.5) << "kd value out of range";
    EXPECT_LE(kd1, 4.5) << "kd value out of range";
    EXPECT_GE(kd2, -4.5) << "kd value out of range";
    EXPECT_LE(kd2, 4.5) << "kd value out of range";
    EXPECT_GE(kd3, -4.5) << "kd value out of range";
    EXPECT_LE(kd3, 4.5) << "kd value out of range";

    std::cout << "kdFromLogP test results:" << std::endl;
    std::cout << "  LogP=0.0 -> kd=" << kd1 << std::endl;
    std::cout << "  LogP=1.0 -> kd=" << kd2 << std::endl;
    std::cout << "  LogP=-1.0 -> kd=" << kd3 << std::endl;
}
