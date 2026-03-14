#include "pharmacophore.h"
#include <GraphMol/Conformer.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "rdkit_utils.h"
#include "surface.h"
#include "rdkit_base_features.h"

PharmacophoreProperty::PharmacophoreProperty(PharmaType type) : type_(type) {
}

std::string PharmacophoreProperty::key() const {
    switch (type_) {
    case PharmaType::ACCEPTOR:
        return "pharma_acceptor";
    case PharmaType::AROMATIC:
        return "pharma_aromatic";
    case PharmaType::DONOR:
        return "pharma_donor";
    case PharmaType::HYDROPHOBE:
        return "pharma_hydrophobe";
    case PharmaType::LUMPED_HYDROPHOBE:
        return "pharma_lumped_hydrophobe";
    case PharmaType::NEG_IONIZABLE:
        return "pharma_neg_ionizable";
    case PharmaType::POS_IONIZABLE:
        return "pharma_pos_ionizable";
    case PharmaType::ZN_BINDER:
        return "pharma_zn_binder";
    }
    return "pharma_unknown";
}

std::string PharmacophoreProperty::getFamilyName() const {
    switch (type_) {
    case PharmaType::ACCEPTOR:
        return "Acceptor";
    case PharmaType::AROMATIC:
        return "Aromatic";
    case PharmaType::DONOR:
        return "Donor";
    case PharmaType::HYDROPHOBE:
        return "Hydrophobe";
    case PharmaType::LUMPED_HYDROPHOBE:
        return "LumpedHydrophobe";
    case PharmaType::NEG_IONIZABLE:
        return "NegIonizable";
    case PharmaType::POS_IONIZABLE:
        return "PosIonizable";
    case PharmaType::ZN_BINDER:
        return "ZnBinder";
    }
    return "";
}

std::vector<PharmacophoreProperty::Feature> PharmacophoreProperty::extractFeatures(
    std::shared_ptr<RDKit::ROMol> mol) const {
    std::vector<Feature> features;
    if (mol->getNumConformers() == 0)
        return features;

    // Use RDKit's built-in feature factory with our embedded definitions
    std::unique_ptr<RDKit::MolChemicalFeatureFactory> factory(
        RDKit::buildFeatureFactory(RDKIT_BASE_FEATURES_FDEF));
    
    if (!factory) {
        std::cerr << "Failed to build RDKit feature factory" << std::endl;
        return features;
    }

    std::string target_family = getFamilyName();
    auto feats = factory->getFeaturesForMol(*mol);
    auto& conf = mol->getConformer();
    
    for (const auto& feat : feats) {
        if (feat->getFamily() == target_family) {
            auto pos = feat->getPos();
            Eigen::Vector3d center(pos.x, pos.y, pos.z);
            
            bool has_direction = false;
            Eigen::Vector3d direction(0, 0, 0);
            std::vector<Eigen::Vector3d> additional_dirs;
            
            // For aromatic rings, try to compute the normal
            if (type_ == PharmaType::AROMATIC) {
                const auto& atoms = feat->getAtoms();
                if (atoms.size() >= 3) {
                    auto p0 = conf.getAtomPos(atoms[0]->getIdx());
                    auto p1 = conf.getAtomPos(atoms[1]->getIdx());
                    auto p2 = conf.getAtomPos(atoms[2]->getIdx());
                    
                    Eigen::Vector3d v1(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                    Eigen::Vector3d v2(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
                    direction = v1.cross(v2).normalized();
                    has_direction = true;
                }
            } else if (type_ == PharmaType::DONOR || type_ == PharmaType::ACCEPTOR) {
                const auto& atoms = feat->getAtoms();
                if (!atoms.empty()) {
                    int idx = atoms[0]->getIdx();
                    int heavy_neighbors = 0;
                    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
                        if (i != static_cast<unsigned int>(idx)) {
                            if (mol->getBondBetweenAtoms(idx, i)) {
                                if (mol->getAtomWithIdx(i)->getSymbol() != "H") {
                                    heavy_neighbors++;
                                }
                            }
                        }
                    }
                    
                    std::pair<std::vector<std::pair<RDKitUtils::Vector3D, RDKitUtils::Vector3D>>, std::string> vects;
                    std::vector<int> featAtoms;
                    for (const auto* a : atoms) {
                        featAtoms.push_back(a->getIdx());
                    }
                    
                    if (type_ == PharmaType::DONOR) {
                        if (heavy_neighbors <= 1) {
                            vects = RDKitUtils::getDonor1FeatVects(conf, featAtoms, 1.0);
                        } else if (heavy_neighbors == 2) {
                            vects = RDKitUtils::getDonor2FeatVects(conf, featAtoms, 1.0);
                        } else {
                            vects = RDKitUtils::getDonor3FeatVects(conf, featAtoms, 1.0);
                        }
                    } else {
                        if (heavy_neighbors <= 1) {
                            vects = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);
                        } else if (heavy_neighbors == 2) {
                            vects = RDKitUtils::getAcceptor2FeatVects(conf, featAtoms, 1.0);
                        } else {
                            vects = RDKitUtils::getAcceptor3FeatVects(conf, featAtoms, 1.0);
                        }
                    }
                    
                    if (!vects.first.empty()) {
                        for (size_t i = 0; i < vects.first.size(); ++i) {
                            Eigen::Vector3d d = (vects.first[i].second.toEigen() - vects.first[i].first.toEigen());
                            if (d.norm() > 1e-5) {
                                if (!has_direction) {
                                    direction = d.normalized();
                                    has_direction = true;
                                } else {
                                    additional_dirs.push_back(d.normalized());
                                }
                            }
                        }
                    }
                }
            }
            
            features.push_back({center, has_direction, direction});
            for (const auto& d : additional_dirs) {
                features.push_back({center, true, d});
            }
        }
    }
    
    return features;
}

void PharmacophoreProperty::compute(const Surface& surface,
                                    std::unordered_map<std::string, std::any>& cache) {
    auto mol_ptr = surface.molecule.get_mol();
    if (!mol_ptr)
        return;

    std::shared_ptr<RDKit::ROMol> mol_copy(new RDKit::ROMol(*mol_ptr));
    std::vector<Feature> features = extractFeatures(mol_copy);

    std::vector<double> values;
    values.reserve(surface.vertices.size());

    double max_dist = 3.0; // Angstroms

    for (const auto& [idx, vertex_data] : surface.vertices) {
        const auto& coords = vertex_data.second;
        Eigen::Vector3d vertex_pos(coords[0], coords[1], coords[2]);

        Eigen::Vector3d surface_normal(0, 0, 0);
        if (coords.size() >= 6) {
            surface_normal = Eigen::Vector3d(coords[3], coords[4], coords[5]).normalized();
        }

        double val = 0.0;
        for (const auto& feat : features) {
            Eigen::Vector3d diff = vertex_pos - feat.position;
            double dist = diff.norm();

            if (dist < max_dist) {
                double feature_val = std::exp(-dist);

                if (feat.has_direction && surface_normal.norm() > 0.1) {
                    Eigen::Vector3d dir_to_surface = diff.normalized();
                    double dot = feat.direction.dot(dir_to_surface);
                    
                    if (type_ == PharmaType::AROMATIC) {
                        feature_val *= std::abs(dot);
                    } else if (type_ == PharmaType::DONOR || type_ == PharmaType::ACCEPTOR) {
                        feature_val *= std::max(0.0, dot);
                    } else {
                        feature_val *= std::abs(dot);
                    }
                }
                
                val = std::max(val, feature_val);
            }
        }
        values.push_back(val);
    }

    std::string key_name = key();
    cache[key_name + "_values"] = values;
    cache[key_name + "_min"] = 0.0;
    cache[key_name + "_max"] = 1.0;
}
