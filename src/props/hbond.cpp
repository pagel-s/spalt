#include "hbond.h"
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "property_params.h"
#include "surface.h"

HBondProperty::HBondProperty() {
}

HBondProperty::~HBondProperty() {
}

void HBondProperty::compute(const Surface& surface,
                            std::unordered_map<std::string, std::any>& cache) {
    try {
        // Get the molecule from the surface
        const Molecule& molecule = surface.molecule;

        // Get the RDKit molecule
        auto mol = molecule.get_mol();
        if (!mol) {
            throw std::runtime_error("Cannot compute hydrogen bond potential: molecule is null");
        }

        // Create a copy for analysis
        auto mol_copy = std::make_shared<RDKit::ROMol>(*mol);

        // Get donor hydrogens and acceptors
        std::set<int> donor_hydrogens = getDonorHydrogens(mol_copy);
        std::set<int> acceptors = getAcceptors(mol_copy);

#ifdef DEBUG
        std::cout << "Found " << donor_hydrogens.size() << " donor hydrogens and "
                  << acceptors.size() << " acceptor atoms" << std::endl;
#endif

        // Extract atom names from surface vertices and compute hydrogen bond potentials
        std::vector<double> hbond_values;
        hbond_values.reserve(surface.vertices.size());

        for (const auto& vertex_pair : surface.vertices) {
            const std::string& atom_name = vertex_pair.second.first;
            const auto& coords = vertex_pair.second.second;

            // Parse atom name according to Python format: "ElementIndex" (name_idx=False)
            std::string atom_key = atom_name;

            // Find the corresponding atom index
            int atom_idx = -1;
            for (unsigned int i = 0; i < mol_copy->getNumAtoms(); ++i) {
                auto atom = mol_copy->getAtomWithIdx(i);
                std::string expected_name = atom->getSymbol() + std::to_string(i);
                if (expected_name == atom_key) {
                    atom_idx = static_cast<int>(i);
                    break;
                }
            }

            if (atom_idx == -1) {
                // If not found, use default value
#ifdef DEBUG
                std::cout << "WARNING: Could not map vertex atom name '" << atom_key
                          << "' to any atom" << std::endl;
#endif
                hbond_values.push_back(0.0);
                continue;
            }

            // Create surface point vector
            Eigen::Vector3d surface_point(coords[0], coords[1], coords[2]);

            // Compute hydrogen bond potential using the core helper function
            double hbond_potential = computeChargeHelperMol(mol_copy, atom_idx, donor_hydrogens,
                                                            acceptors, surface_point);
            hbond_values.push_back(hbond_potential);
        }

        // Store the results in cache
        cache["hbond_values"] = hbond_values;

        // Compute statistics for display
        if (!hbond_values.empty()) {
            double avg_hbond = std::accumulate(hbond_values.begin(), hbond_values.end(), 0.0) /
                               hbond_values.size();

            cache["hbond_avg"] = avg_hbond;
#ifdef DEBUG
            double min_hbond = *std::min_element(hbond_values.begin(), hbond_values.end());
            double max_hbond = *std::max_element(hbond_values.begin(), hbond_values.end());
            std::cout << "Hydrogen bond potential computation successful: " << hbond_values.size()
                      << " values" << std::endl;
            std::cout << "  Actual min hbond potential: " << min_hbond << std::endl;
            std::cout << "  Actual max hbond potential: " << max_hbond << std::endl;
            std::cout << "  Average hbond potential: " << avg_hbond << std::endl;
#endif
        }

        // Store standardized min/max values for consistent color mapping across molecules
        PropertyParams hbond_params = PropertyParamsRegistry::getParams("hbond");
        cache["hbond_min"] = hbond_params.min_value;
        cache["hbond_max"] = hbond_params.max_value;

    } catch (const std::exception& e) {
        std::cerr << "Error computing hydrogen bond potential: " << e.what() << std::endl;
        throw;
    }
}

std::set<int> HBondProperty::getDonorHydrogens(std::shared_ptr<RDKit::ROMol> mol) const {
    std::set<int> donor_hydrogens;

    try {
        std::set<int> donors = getDonors(mol);

        // Find hydrogens attached to donor atoms
        for (int donor_idx : donors) {
            // Check all atoms to find hydrogens bonded to this donor
            for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
                auto other_atom = mol->getAtomWithIdx(j);
                if (other_atom->getSymbol() == "H") {
                    auto bond = mol->getBondBetweenAtoms(donor_idx, j);
                    if (bond) {
                        donor_hydrogens.insert(j);
                    }
                }
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error getting donor hydrogens: " << e.what() << std::endl;
    }

    return donor_hydrogens;
}

std::set<int> HBondProperty::getAcceptors(std::shared_ptr<RDKit::ROMol> mol) const {
    std::set<int> acceptors;

    try {
        // Define SMARTS patterns for acceptors based on the fdef file
        // These match the SingleAtomAcceptor feature definition
        std::vector<std::string> acceptor_smarts = {
            // Hydroxyl acceptors: [O;H1;v2]
            "[O;H1;v2]",

            // ChalcAcceptor patterns:
            "[O;H0;v2;!$(O=N-*)]", "[O;-;!$(*-N=O)]", "[o;+0]",

            // NAcceptor patterns:
            "[n;+0;!X3;!$([n;H1](cc)cc)]", "[N;H0]#[C&v4]", "[N&v3;H0;$(Nc)]",

            // HalogenAcceptor
            "[F;$(F-[#6]);!$(FC[F,Cl,Br,I])]"};

        // Apply each SMARTS pattern using proper RDKit C++ API
        for (const auto& smarts : acceptor_smarts) {
            try {
                auto query = RDKit::SmartsToMol(smarts);
                if (!query) {
                    std::cerr << "Warning: Invalid SMARTS pattern: " << smarts << std::endl;
                    continue;
                }

                // Get all matches at once
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(*mol, *query, matches);
                // std::cout << "Acceptors: Found " << matches.size() << " matches for SMARTS
                // pattern " << smarts << std::endl;
                for (const auto& match : matches) {
                    if (!match.empty()) {
                        acceptors.insert(
                            match[0].second);  // Add atom index (second element of pair)
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Warning: Error processing SMARTS pattern " << smarts << ": "
                          << e.what() << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error getting acceptors: " << e.what() << std::endl;
    }

    return acceptors;
}

std::set<int> HBondProperty::getDonors(std::shared_ptr<RDKit::ROMol> mol) const {
    std::set<int> donors;

    try {
        // Define SMARTS patterns for donors based on the fdef file
        // These match the SingleAtomDonor feature definition
        std::vector<std::string> donor_smarts = {
            // NDonor patterns
            "[N&!H0&v3,N&!H0&+1&v4,n&H1&+0]", "[$([Nv3](-C)(-C)-C)]", "[$(n[n;H1]),$(nc[n;H1])]",
            // ChalcDonor patterns
            "[O,S;H1;+0]"};

        // Apply each SMARTS pattern using proper RDKit C++ API
        for (const auto& smarts : donor_smarts) {
            try {
                auto query = RDKit::SmartsToMol(smarts);
                if (!query) {
                    std::cerr << "Warning: Invalid SMARTS pattern: " << smarts << std::endl;
                    continue;
                }

                // Get all matches at once
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(*mol, *query, matches);

                // std::cout << "Donors: Found " << matches.size() << " matches for SMARTS pattern "
                // << smarts << std::endl;
                for (const auto& match : matches) {
                    if (!match.empty()) {
                        donors.insert(match[0].second);  // Add atom index (second element of pair)
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Warning: Error processing SMARTS pattern " << smarts << ": "
                          << e.what() << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error getting donors: " << e.what() << std::endl;
    }

    return donors;
}

int HBondProperty::getAcceptorType(const RDKit::Atom& acceptor_atom) const {
    // This matches Python: return sum([a.GetSymbol() != 'H' for a in acceptor_atom.GetNeighbors()])
    // Acceptor of type 1: acceptor with one heavy atom neighbor.
    // Acceptor of type 2: acceptor with two adjacent heavy atoms.
    // Acceptor of type 3: acceptor with three adjacent heavy atoms.

    int heavy_neighbors = 0;
    // Count heavy atom neighbors (non-hydrogen) - simplified approach since we don't have
    // GetNeighbors()
    for (unsigned int i = 0; i < acceptor_atom.getOwningMol().getNumAtoms(); ++i) {
        if (i != acceptor_atom.getIdx()) {
            auto bond = acceptor_atom.getOwningMol().getBondBetweenAtoms(acceptor_atom.getIdx(), i);
            if (bond) {
                auto neighbor_atom = acceptor_atom.getOwningMol().getAtomWithIdx(i);
                if (neighbor_atom->getSymbol() != "H") {
                    heavy_neighbors++;
                }
            }
        }
    }
    return heavy_neighbors;
}

double HBondProperty::computeAngleDeviation(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                            const Eigen::Vector3d& c, double ideal_angle) const {
    // Compute vectors from b to a and b to c
    Eigen::Vector3d ba = a - b;
    Eigen::Vector3d bc = c - b;

    // Normalize vectors
    ba.normalize();
    bc.normalize();

    // Compute angle between vectors
    double cos_angle = ba.dot(bc);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));  // Clamp to avoid numerical errors
    double actual_angle = std::acos(cos_angle);

    // Return absolute deviation from ideal angle
    return std::abs(actual_angle - ideal_angle);
}

double HBondProperty::computeAnglePenalty(double angle_deviation) const {
    // Standard deviation: hbond_std_dev
    double penalty = 1.0 - std::pow(angle_deviation / HBOND_STD_DEV, 2);
    return std::max(0.0, penalty);
}

double HBondProperty::computeChargeHelperMol(std::shared_ptr<RDKit::ROMol> mol,
                                             unsigned int atom_idx,
                                             const std::set<int>& donor_hydrogens,
                                             const std::set<int>& acceptors,
                                             const Eigen::Vector3d& surface_point) const {
    // Check if it is a polar hydrogen (donor)
    if (donor_hydrogens.find(atom_idx) != donor_hydrogens.end()) {
        int donor_idx = -1;
        for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
            if (j != atom_idx) {
                auto bond = mol->getBondBetweenAtoms(atom_idx, j);
                if (bond) {
                    auto other_atom = mol->getAtomWithIdx(j);
                    if (other_atom->getSymbol() != "H") {
                        donor_idx = j;
                        break;
                    }
                }
            }
        }

        if (donor_idx != -1) {
            // Get coordinates - matches Python: a =
            // rdmol.GetConformer().GetAtomPosition(donor.GetIdx())
            auto conf = mol->getConformer(0);
            Eigen::Vector3d donor_pos(conf.getAtomPos(donor_idx).x, conf.getAtomPos(donor_idx).y,
                                      conf.getAtomPos(donor_idx).z);
            // b = rdmol.GetConformer().GetAtomPosition(atom_idx)
            Eigen::Vector3d hydrogen_pos(conf.getAtomPos(atom_idx).x, conf.getAtomPos(atom_idx).y,
                                         conf.getAtomPos(atom_idx).z);

            // Donor-H is always 180.0 degrees (π) - matches Python: angle_deviation =
            // computeAngleDeviation(a, b, v, np.pi)
            double angle_deviation =
                computeAngleDeviation(donor_pos, hydrogen_pos, surface_point, M_PI);
            double angle_penalty = computeAnglePenalty(angle_deviation);
            return 1.0 * angle_penalty;  // matches Python: return 1.0 * angle_penalty
        }
    }

    // Check if it is an acceptor oxygen or nitrogen - matches Python: elif atom_idx in acceptors:
    else if (acceptors.find(atom_idx) != acceptors.end()) {
        auto vects_pair = getAcceptorFeatVects(mol, atom_idx);
        auto vects = vects_pair.first;
        std::string vec_type = vects_pair.second;

        // Note: the ideal direction is often ambiguous, e.g. one of two options if two lone
        // electron pairs are present or anywhere on a cone. We compute the smallest possible angle
        // deviation in these cases.
        double angle_deviation;
        if (vec_type == "cone") {
            // the ideal direction is anywhere on a cone at 60 degrees from this vector
            // angle_deviation = computeAngleDeviation(ep, sp, v, np.pi / 3)
            auto sp_ep = vects[0];
            Eigen::Vector3d sp = sp_ep.first;
            Eigen::Vector3d ep = sp_ep.second;
            angle_deviation = computeAngleDeviation(ep, sp, surface_point, M_PI / 3.0);
        } else {
            // consider all options and select minimum
            // angle_deviation = min([computeAngleDeviation(ep, sp, v, 0) for sp, ep in vects])
            double min_angle_deviation = std::numeric_limits<double>::max();
            for (const auto& sp_ep : vects) {
                Eigen::Vector3d sp = sp_ep.first;
                Eigen::Vector3d ep = sp_ep.second;
                double current_angle_deviation = computeAngleDeviation(ep, sp, surface_point, 0.0);
                min_angle_deviation = std::min(min_angle_deviation, current_angle_deviation);
            }
            angle_deviation = min_angle_deviation;
        }

        // angle_penalty = computeAnglePenalty(angle_deviation)
        double angle_penalty = computeAnglePenalty(angle_deviation);

        // There's no need to compute the plane penalty separately, out-of-plane directions are
        // already penalized plane_penalty = 1.0
        double plane_penalty = 1.0;

        // return -1.0 * angle_penalty * plane_penalty
        return -1.0 * angle_penalty * plane_penalty;
    }

    return 0.0;
}

double HBondProperty::computeHBondPotential(std::shared_ptr<RDKit::ROMol> mol, int atom_idx,
                                            const std::set<int>& donor_hydrogens,
                                            const std::set<int>& acceptors,
                                            const Eigen::Vector3d& surface_point) const {
    // Legacy function - now calls the correct computeChargeHelperMol
    return computeChargeHelperMol(mol, atom_idx, donor_hydrogens, acceptors, surface_point);
}

std::pair<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>, std::string>
HBondProperty::getAcceptorFeatVects(std::shared_ptr<RDKit::ROMol> mol, int atom_idx) const {
    // Use the new RDKit utilities for accurate acceptor direction calculation
    auto conf = mol->getConformer(0);
    std::vector<int> featAtoms = {atom_idx};

    // Get acceptor type to determine which function to call
    // TODO: this will have to be refined to find acceptors properly!
    auto acceptor_atom = mol->getAtomWithIdx(atom_idx);
    int acceptor_type = getAcceptorType(*acceptor_atom);

    // Call the appropriate RDKit utility function based on acceptor type
    std::pair<std::vector<std::pair<RDKitUtils::Vector3D, RDKitUtils::Vector3D>>, std::string>
        result;

    if (acceptor_type == 1) {
        result = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);
    } else if (acceptor_type == 2) {
        result = RDKitUtils::getAcceptor2FeatVects(conf, featAtoms, 1.0);
    } else if (acceptor_type == 3) {
        result = RDKitUtils::getAcceptor3FeatVects(conf, featAtoms, 1.0);
    } else {
        // For acceptor types > 3, return empty result
        result = std::make_pair(
            std::vector<std::pair<RDKitUtils::Vector3D, RDKitUtils::Vector3D>>(), "none");
    }

    // Convert from RDKitUtils::Vector3D to Eigen::Vector3d
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> eigen_result;
    for (const auto& pair : result.first) {
        Eigen::Vector3d start(pair.first.x, pair.first.y, pair.first.z);
        Eigen::Vector3d end(pair.second.x, pair.second.y, pair.second.z);
        eigen_result.push_back(std::make_pair(start, end));
    }

    return std::make_pair(eigen_result, result.second);
}
