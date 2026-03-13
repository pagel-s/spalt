#include <gtest/gtest.h>
#include "molecule.h"

class PharmacophoreTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Benzene (aromatic)
        mol_benzene = std::make_unique<Molecule>("C1=CC=CC=C1");
        mol_benzene->calculate_conformers(1);
        // Methylamine (donor)
        mol_amine = std::make_unique<Molecule>("CN");
        mol_amine->calculate_conformers(1);
        // Acetone (acceptor)
        mol_ketone = std::make_unique<Molecule>("CC(=O)C");
        mol_ketone->calculate_conformers(1);
    }

    std::unique_ptr<Molecule> mol_benzene;
    std::unique_ptr<Molecule> mol_amine;
    std::unique_ptr<Molecule> mol_ketone;
};

TEST_F(PharmacophoreTest, ComputePharmacophoreAromatic) {
    mol_benzene->compute("pharma_aromatic", 0);
    auto vals = mol_benzene->getProperty<std::vector<double>>("pharma_aromatic",
                                                              "pharma_aromatic_values", 0);
    EXPECT_FALSE(vals.empty());

    // Benzene should have some positive aromatic values
    bool has_aromatic = false;
    for (double v : vals) {
        if (v > 0.0)
            has_aromatic = true;
    }
    EXPECT_TRUE(has_aromatic);
}

TEST_F(PharmacophoreTest, ComputePharmacophoreDonor) {
    mol_amine->compute("pharma_donor", 0);
    auto vals =
        mol_amine->getProperty<std::vector<double>>("pharma_donor", "pharma_donor_values", 0);
    EXPECT_FALSE(vals.empty());
}

TEST_F(PharmacophoreTest, ComputePharmacophoreAcceptor) {
    mol_ketone->compute("pharma_acceptor", 0);
    auto vals = mol_ketone->getProperty<std::vector<double>>("pharma_acceptor",
                                                             "pharma_acceptor_values", 0);
    EXPECT_FALSE(vals.empty());
}
