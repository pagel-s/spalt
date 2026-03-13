#include "hydrophobicity.h"
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <algorithm>
#include <boost/make_shared.hpp>
#include <cmath>
#include <iostream>
#include "property_params.h"
#include "surface.h"

HydrophobicityProperty::HydrophobicityProperty() {
}

HydrophobicityProperty::~HydrophobicityProperty() {
}

void HydrophobicityProperty::compute(const Surface& surface,
                                     std::unordered_map<std::string, std::any>& cache) {
    try {
        // Get the molecule from the surface
        const Molecule& molecule = surface.molecule;

        // Get the RDKit molecule
        auto mol = molecule.get_mol();
        if (!mol) {
            throw std::runtime_error("Cannot compute hydrophobicity: molecule is null");
        }

        // Create a copy for fragmentation (we'll modify it)
        auto mol_copy = std::make_shared<RDKit::ROMol>(*mol);

        // Assign hydrophobicity values to atoms
        std::map<std::string, double> atom_kd = molAssignKd(mol_copy);

        // Extract atom names from surface vertices and compute hydrophobicity
        std::vector<double> hydrophobicity_values;
        hydrophobicity_values.reserve(surface.vertices.size());

        for (const auto& vertex_pair : surface.vertices) {
            const std::string& atom_name = vertex_pair.second.first;

            // Extract the atom identifier from the name (format: "ElementIndex")
            // The name format should match what's generated in generateXYZRFile
            std::string atom_key = atom_name;

            // Try to find the hydrophobicity value for this atom
            auto it = atom_kd.find(atom_key);
            if (it != atom_kd.end()) {
                hydrophobicity_values.push_back(it->second);
            } else {
                // If not found, use a default value (neutral hydrophobicity)
                hydrophobicity_values.push_back(0.0);
#ifdef DEBUG
                std::cout << "Warning: No hydrophobicity value found for atom " << atom_key
                          << ", using default 0.0" << std::endl;
#endif
            }
        }

        // Store the results in cache
        cache["hydrophobicity_values"] = hydrophobicity_values;

        // Compute statistics for display
        if (!hydrophobicity_values.empty()) {
            double avg_hyd =
                std::accumulate(hydrophobicity_values.begin(), hydrophobicity_values.end(), 0.0) /
                hydrophobicity_values.size();

            cache["hydrophobicity_avg"] = avg_hyd;
#ifdef DEBUG
            double min_hyd =
                *std::min_element(hydrophobicity_values.begin(), hydrophobicity_values.end());
            double max_hyd =
                *std::max_element(hydrophobicity_values.begin(), hydrophobicity_values.end());
            std::cout << "Hydrophobicity computation successful: " << hydrophobicity_values.size()
                      << " values" << std::endl;
            std::cout << "  Actual min hydrophobicity: " << min_hyd << std::endl;
            std::cout << "  Actual max hydrophobicity: " << max_hyd << std::endl;
            std::cout << "  Average hydrophobicity: " << avg_hyd << std::endl;
#endif
        }

        // Store standardized min/max values for consistent color mapping across molecules
        PropertyParams hydrophobicity_params = PropertyParamsRegistry::getParams("hydrophobicity");
        cache["hydrophobicity_min"] = hydrophobicity_params.min_value;
        cache["hydrophobicity_max"] = hydrophobicity_params.max_value;

    } catch (const std::exception& e) {
        std::cerr << "Error computing hydrophobicity: " << e.what() << std::endl;
        throw;
    }
}

double HydrophobicityProperty::kdFromLogP(double logp, double kd_min, double kd_max) const {
    // Convert LogP to hydrophobicity using the formula from the Python reference
    double kd = -6.2786 + std::exp(0.4772 * logp + 1.8491);

    // Clip to the specified range
    return std::max(kd_min, std::min(kd_max, kd));
}

std::map<std::string, double> HydrophobicityProperty::molAssignKd(
    std::shared_ptr<RDKit::ROMol> mol) const {
    std::map<std::string, double> atom_kd;

    try {
        // Create atom name mapping (similar to Python version)
        std::map<unsigned int, std::string> atom_name_map;
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            auto atom = mol->getAtomWithIdx(i);
            atom->setAtomMapNum(i + 1);  // 0 is reserved as default value
            std::string atom_name = atom->getSymbol() + std::to_string(i);
            atom_name_map[i + 1] = atom_name;
        }

        // Get fragments using BRICS fragmentation
        std::vector<RDKit::ROMOL_SPTR> fragments = getFragments(mol);

        // Assign hydrophobicity to each atom based on its fragment
        for (const auto& frag : fragments) {
            // Calculate LogP for this fragment
            double frag_logp = RDKit::Descriptors::calcClogP(*frag);
            double frag_kd = kdFromLogP(frag_logp);

            // Assign this kd value to all atoms in the fragment
            for (auto atom : frag->atoms()) {
                unsigned int atom_mapnum = atom->getAtomMapNum();
                if (atom_mapnum == 0) {
                    continue;  // Skip atoms that don't exist in original molecule
                }

                std::string atom_name = atom_name_map[atom_mapnum];
                atom_kd[atom_name] = frag_kd;
            }
        }

        return atom_kd;

    } catch (const std::exception& e) {
        std::cerr << "Error assigning hydrophobicity values: " << e.what() << std::endl;
        throw;
    }
}

std::vector<RDKit::ROMOL_SPTR> HydrophobicityProperty::getFragments(
    std::shared_ptr<RDKit::ROMol> mol) const {
    std::vector<RDKit::ROMOL_SPTR> fragments;

    try {
        // Fragment the molecule on BRICS bonds
        RDKit::ROMol* fragmented_mol = RDKit::MolFragmenter::fragmentOnBRICSBonds(*mol);

        if (!fragmented_mol) {
            throw std::runtime_error("Failed to fragment molecule on BRICS bonds");
        }

        // Get molecular fragments
        std::vector<RDKit::ROMOL_SPTR> mol_frags =
            RDKit::MolOps::getMolFrags(*fragmented_mol, true);

        // Find exit atoms (dummy atoms with symbol '*') to be removed
        for (auto& frag : mol_frags) {
            std::vector<unsigned int> exits;
            for (auto atom : frag->atoms()) {
                if (atom->getSymbol() == "*") {
                    exits.push_back(atom->getIdx());
                }
            }

            // Remove exit atoms in reverse order to maintain indices
            std::sort(exits.rbegin(), exits.rend());
            RDKit::RWMol editable_frag(*frag);
            for (unsigned int idx : exits) {
                editable_frag.removeAtom(idx);
            }

            // Add the cleaned fragment
            fragments.push_back(boost::make_shared<RDKit::ROMol>(editable_frag));
        }

        // Clean up
        delete fragmented_mol;

        return fragments;

    } catch (const std::exception& e) {
        std::cerr << "Error fragmenting molecule: " << e.what() << std::endl;
        throw;
    }
}
