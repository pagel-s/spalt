#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include "molecule.h"
#include "surface.h"
#include "props/hbond.h"

TEST(VertexMappingTest, TestVertexToAtomMapping) {
    // Load test molecule
    Molecule mol("data/Mol1.sdf");
    
    // Get RDKit molecule and convert to shared_ptr
    auto rdmol_ptr = mol.get_mol();
    ASSERT_TRUE(rdmol_ptr) << "Failed to get RDKit molecule";
    std::shared_ptr<RDKit::ROMol> rdmol(rdmol_ptr.release());
    
    std::cout << "Loaded molecule with " << rdmol->getNumAtoms() << " atoms" << std::endl;
    
    // Generate surface
    Surface surface(mol, 0, 50, 1.2, 3.0, 3.0, "tses", "msms", "full");
    
    std::cout << "Generated surface with " << surface.vertices.size() << " vertices" << std::endl;
    
    HBondProperty hbond_prop;
    auto acceptors = hbond_prop.getAcceptors(rdmol);
    auto donors = hbond_prop.getDonors(rdmol);
    auto donor_hydrogens = hbond_prop.getDonorHydrogens(rdmol);
    
    std::cout << "\n=== DIRECT HBOND DETECTION ===" << std::endl;
    std::cout << "Acceptors: " << acceptors.size() << std::endl;
    for (int acc : acceptors) {
        auto atom = rdmol->getAtomWithIdx(acc);
        std::cout << "  Acceptor " << acc << ": " << atom->getSymbol() << std::endl;
    }
    
    std::cout << "Donors: " << donors.size() << std::endl;
    for (int donor : donors) {
        auto atom = rdmol->getAtomWithIdx(donor);
        std::cout << "  Donor " << donor << ": " << atom->getSymbol() << std::endl;
    }
    
    std::cout << "Donor hydrogens: " << donor_hydrogens.size() << std::endl;
    for (int h : donor_hydrogens) {
        auto atom = rdmol->getAtomWithIdx(h);
        std::cout << "  Donor H " << h << ": " << atom->getSymbol() << std::endl;
    }
    
    // Test the surface vertex mapping
    std::cout << "\n=== SURFACE VERTEX MAPPING ===" << std::endl;
    
    // Check first 10 vertices and their atom names
    int count = 0;
    for (const auto& vertex_pair : surface.vertices) {
        if (count >= 10) break;
        
        const std::string& atom_name = vertex_pair.second.first;
        std::cout << "Vertex " << vertex_pair.first << " -> Atom name: '" << atom_name << "'" << std::endl;
        
        // Try to extract atom index from name
        try {
            // Assuming format like "ElementIndex" or similar
            size_t last_digit = atom_name.find_last_not_of("0123456789");
            if (last_digit != std::string::npos && last_digit + 1 < atom_name.length()) {
                std::string index_str = atom_name.substr(last_digit + 1);
                int atom_idx = std::stoi(index_str);
                
                if (atom_idx < static_cast<int>(rdmol->getNumAtoms())) {
                    auto atom = rdmol->getAtomWithIdx(atom_idx);
                    std::cout << "  -> Mapped to atom " << atom_idx << ": " << atom->getSymbol();
                    
                    // Check if this atom is a donor/acceptor
                    if (acceptors.count(atom_idx)) {
                        std::cout << " (ACCEPTOR)";
                    }
                    if (donors.count(atom_idx)) {
                        std::cout << " (DONOR)";
                    }
                    std::cout << std::endl;
                } else {
                    std::cout << "  -> INVALID atom index: " << atom_idx << std::endl;
                }
            } else {
                std::cout << "  -> Could not extract atom index from name" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  -> Error parsing atom name: " << e.what() << std::endl;
        }
        
        count++;
    }
    
    // Test the actual hydrogen bond property computation on surface
    std::cout << "\n=== SURFACE HBOND COMPUTATION ===" << std::endl;
    std::unordered_map<std::string, std::any> cache;
    hbond_prop.compute(surface, cache);
    
    // Check if we got any hbond values
    auto hbond_it = cache.find("hbond");
    ASSERT_NE(hbond_it, cache.end()) << "No hbond values found in cache!";
    
    try {
        auto hbond_values = std::any_cast<std::vector<double>>(hbond_it->second);
        std::cout << "Computed " << hbond_values.size() << " hbond values" << std::endl;
        
        // Count non-zero values
        int non_zero_count = 0;
        double max_val = 0.0;
        for (double val : hbond_values) {
            if (std::abs(val) > 1e-6) {
                non_zero_count++;
                max_val = std::max(max_val, std::abs(val));
            }
        }
        std::cout << "Non-zero hbond values: " << non_zero_count << std::endl;
        std::cout << "Max absolute hbond value: " << max_val << std::endl;
        
        // This is the key test - we should have some non-zero values if mapping is working
        EXPECT_GT(non_zero_count, 0) << "No non-zero hydrogen bond values found - vertex-to-atom mapping may be broken";
        
    } catch (const std::exception& e) {
        FAIL() << "Error accessing hbond values: " << e.what();
    }
}
