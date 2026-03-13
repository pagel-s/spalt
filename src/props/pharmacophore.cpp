#include "pharmacophore.h"
#include <GraphMol/Conformer.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RingInfo.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "hbond.h"
#include "rdkit_utils.h"
#include "surface.h"

PharmacophoreProperty::PharmacophoreProperty(PharmaType type) : type_(type) {
}

std::string PharmacophoreProperty::key() const {
    switch (type_) {
    case PharmaType::AROMATIC:
        return "pharma_aromatic";
    case PharmaType::POSITIVE:
        return "pharma_pos";
    case PharmaType::NEGATIVE:
        return "pharma_neg";
    case PharmaType::DONOR:
        return "pharma_donor";
    case PharmaType::ACCEPTOR:
        return "pharma_acceptor";
    }
    return "pharma_unknown";
}

std::vector<PharmacophoreProperty::Feature> PharmacophoreProperty::extractFeatures(
    std::shared_ptr<RDKit::ROMol> mol) const {
    std::vector<Feature> features;
    if (mol->getNumConformers() == 0)
        return features;
    auto& conf = mol->getConformer();

    if (type_ == PharmaType::AROMATIC) {
        auto ringInfo = mol->getRingInfo();
        if (ringInfo && ringInfo->isInitialized()) {
            for (const auto& ring : ringInfo->atomRings()) {
                bool isAromatic = true;
                for (int idx : ring) {
                    if (!mol->getAtomWithIdx(idx)->getIsAromatic()) {
                        isAromatic = false;
                        break;
                    }
                }
                if (isAromatic && ring.size() >= 3) {
                    Eigen::Vector3d center(0, 0, 0);
                    for (int idx : ring) {
                        auto pos = conf.getAtomPos(idx);
                        center += Eigen::Vector3d(pos.x, pos.y, pos.z);
                    }
                    center /= ring.size();

                    auto pos0 = conf.getAtomPos(ring[0]);
                    auto pos1 = conf.getAtomPos(ring[1]);
                    auto pos2 = conf.getAtomPos(ring[2]);
                    Eigen::Vector3d v1(pos1.x - pos0.x, pos1.y - pos0.y, pos1.z - pos0.z);
                    Eigen::Vector3d v2(pos2.x - pos0.x, pos2.y - pos0.y, pos2.z - pos0.z);
                    Eigen::Vector3d normal = v1.cross(v2).normalized();

                    features.push_back({center, true, normal});
                }
            }
        }
    } else if (type_ == PharmaType::POSITIVE) {
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            auto atom = mol->getAtomWithIdx(i);
            if (atom->getFormalCharge() > 0) {
                auto pos = conf.getAtomPos(i);
                features.push_back(
                    {Eigen::Vector3d(pos.x, pos.y, pos.z), false, Eigen::Vector3d::Zero()});
            }
        }
    } else if (type_ == PharmaType::NEGATIVE) {
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            auto atom = mol->getAtomWithIdx(i);
            if (atom->getFormalCharge() < 0) {
                auto pos = conf.getAtomPos(i);
                features.push_back(
                    {Eigen::Vector3d(pos.x, pos.y, pos.z), false, Eigen::Vector3d::Zero()});
            }
        }
    } else if (type_ == PharmaType::DONOR) {
        HBondProperty hbond;
        std::set<int> donors = hbond.getDonors(mol);
        for (int idx : donors) {
            auto pos = conf.getAtomPos(idx);
            Eigen::Vector3d center(pos.x, pos.y, pos.z);

            bool found_dir = false;
            Eigen::Vector3d dir(0, 0, 0);

            for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
                auto other_atom = mol->getAtomWithIdx(j);
                if (other_atom->getSymbol() == "H") {
                    if (mol->getBondBetweenAtoms(idx, j)) {
                        auto h_pos = conf.getAtomPos(j);
                        dir += Eigen::Vector3d(h_pos.x - pos.x, h_pos.y - pos.y, h_pos.z - pos.z)
                                   .normalized();
                        found_dir = true;
                    }
                }
            }
            if (found_dir) {
                features.push_back({center, true, dir.normalized()});
            } else {
                features.push_back({center, false, Eigen::Vector3d::Zero()});
            }
        }
    } else if (type_ == PharmaType::ACCEPTOR) {
        HBondProperty hbond;
        std::set<int> acceptors = hbond.getAcceptors(mol);
        for (int idx : acceptors) {
            auto pos = conf.getAtomPos(idx);
            Eigen::Vector3d center(pos.x, pos.y, pos.z);

            std::vector<int> featAtoms = {idx};
            int type = hbond.getAcceptorType(*mol->getAtomWithIdx(idx));

            std::pair<std::vector<std::pair<RDKitUtils::Vector3D, RDKitUtils::Vector3D>>,
                      std::string>
                vects;
            if (type == 1) {
                vects = RDKitUtils::getAcceptor1FeatVects(conf, featAtoms, 1.0);
            } else if (type == 2) {
                vects = RDKitUtils::getAcceptor2FeatVects(conf, featAtoms, 1.0);
            } else {
                vects = RDKitUtils::getAcceptor3FeatVects(conf, featAtoms, 1.0);
            }

            if (!vects.first.empty()) {
                Eigen::Vector3d dir(0, 0, 0);
                for (const auto& pair : vects.first) {
                    dir += (pair.second.toEigen() - pair.first.toEigen()).normalized();
                }
                if (dir.norm() > 1e-5) {
                    features.push_back({center, true, dir.normalized()});
                } else {
                    features.push_back({center, false, Eigen::Vector3d::Zero()});
                }
            } else {
                features.push_back({center, false, Eigen::Vector3d::Zero()});
            }
        }
    }

    return features;
}

void PharmacophoreProperty::compute(const Surface& surface,
                                    std::unordered_map<std::string, std::any>& cache) {
    auto mol = surface.molecule.get_mol();
    if (!mol) {
        throw std::runtime_error("Cannot compute pharmacophore potential: molecule is null");
    }
    auto mol_copy = std::make_shared<RDKit::ROMol>(*mol);
    std::vector<Feature> features = extractFeatures(mol_copy);

    std::vector<double> values;
    values.reserve(surface.vertices.size());

    for (const auto& vertex_pair : surface.vertices) {
        Eigen::Vector3d p(vertex_pair.second.second[0], vertex_pair.second.second[1],
                          vertex_pair.second.second[2]);

        double max_val = 0.0;
        for (const auto& feat : features) {
            double d = (p - feat.position).norm();
            // Faster decay so features don't bleed too much.
            double dist_penalty = std::exp(-0.5 * d * d);

            double val = dist_penalty;
            if (feat.has_direction) {
                Eigen::Vector3d v_to_p = (p - feat.position).normalized();
                if (type_ == PharmaType::AROMATIC) {
                    // Aromatic rings have pi clouds above and below
                    double angle_align = std::abs(v_to_p.dot(feat.direction));
                    val *= angle_align;
                } else if (type_ == PharmaType::DONOR || type_ == PharmaType::ACCEPTOR) {
                    // Directional preference for H-bonds
                    double angle_align = v_to_p.dot(feat.direction);
                    if (angle_align < 0)
                        angle_align = 0;  // strict forward
                    val *= angle_align;
                }
            }
            max_val = std::max(max_val, val);
        }
        values.push_back(max_val);
    }

    std::string key_name = key();
    cache[key_name + "_values"] = values;
    cache[key_name + "_min"] = 0.0;
    cache[key_name + "_max"] = 1.0;
}
