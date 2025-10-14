#include "rdkit_utils.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RDKitUtils {

// Vector3D implementation
Vector3D Vector3D::operator+(const Vector3D& other) const {
    return Vector3D(x + other.x, y + other.y, z + other.z);
}

Vector3D Vector3D::operator-(const Vector3D& other) const {
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

Vector3D Vector3D::operator*(double scalar) const {
    return Vector3D(x * scalar, y * scalar, z * scalar);
}

Vector3D Vector3D::operator/(double scalar) const {
    if (std::abs(scalar) < 1e-10) {
        throw std::runtime_error("Division by zero");
    }
    return Vector3D(x / scalar, y / scalar, z / scalar);
}

Vector3D Vector3D::operator-() const {
    return Vector3D(-x, -y, -z);
}

double Vector3D::dot(const Vector3D& other) const {
    return x * other.x + y * other.y + z * other.z;
}

Vector3D Vector3D::cross(const Vector3D& other) const {
    return Vector3D(y * other.z - z * other.y, -x * other.z + z * other.x,
                    x * other.y - y * other.x);
}

void Vector3D::normalize() {
    double n = norm();
    if (n > 1e-10) {
        x /= n;
        y /= n;
        z /= n;
    }
}

Vector3D Vector3D::normalized() const {
    Vector3D result = *this;
    result.normalize();
    return result;
}

double Vector3D::norm() const {
    return std::sqrt(x * x + y * y + z * z);
}

double Vector3D::normsq() const {
    return x * x + y * y + z * z;
}

double Vector3D::angle(const Vector3D& other) const {
    double n1 = norm();
    double n2 = other.norm();
    if (n1 < 1e-10 || n2 < 1e-10) {
        return 0.0;
    }

    double c = dot(other) / (n1 * n2);
    // Take care of roundoff errors
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
}

// Cross product function
Vector3D cross(const Vector3D& v1, const Vector3D& v2) {
    return v1.cross(v2);
}

// Find neighbors from adjacency matrix
std::vector<int> findNeighbors(int atomId, const std::vector<std::vector<int>>& adjMat) {
    std::vector<int> nbrs;
    if (atomId >= 0 && atomId < static_cast<int>(adjMat.size())) {
        for (size_t i = 0; i < adjMat[atomId].size(); ++i) {
            if (adjMat[atomId][i] >= 1) {
                nbrs.push_back(static_cast<int>(i));
            }
        }
    }
    return nbrs;
}

// Find average vector
Vector3D findAvgVec(const RDKit::Conformer& conf, const Vector3D& center,
                    const std::vector<const RDKit::Atom*>& nbrs) {
    Vector3D avgVec(0, 0, 0);
    bool first = true;

    for (auto* nbr : nbrs) {
        if (nbr == nullptr)
            continue;

        int nid = nbr->getIdx();
        auto pos = conf.getAtomPos(nid);
        Vector3D pt(pos.x, pos.y, pos.z);
        pt = pt - center;
        pt.normalize();

        if (first) {
            avgVec = pt;
            first = false;
        } else {
            avgVec = avgVec + pt;
        }
    }

    avgVec.normalize();
    return avgVec;
}

// Arbitrary axis rotation
Vector3D arbAxisRotation(double theta, const Vector3D& ax, const Vector3D& pt) {
    double theta_rad = M_PI * theta / 180.0;
    double c = std::cos(theta_rad);
    double s = std::sin(theta_rad);
    double t = 1.0 - c;

    double X = ax.x;
    double Y = ax.y;
    double Z = ax.z;

    // Rotation matrix
    double mat[3][3] = {{t * X * X + c, t * X * Y + s * Z, t * X * Z - s * Y},
                        {t * X * Y - s * Z, t * Y * Y + c, t * Y * Z + s * X},
                        {t * X * Z + s * Y, t * Y * Z - s * X, t * Z * Z + c}};

    // Apply rotation matrix
    double new_x = mat[0][0] * pt.x + mat[0][1] * pt.y + mat[0][2] * pt.z;
    double new_y = mat[1][0] * pt.x + mat[1][1] * pt.y + mat[1][2] * pt.z;
    double new_z = mat[2][0] * pt.x + mat[2][1] * pt.y + mat[2][2] * pt.z;

    return Vector3D(new_x, new_y, new_z);
}

// Get aromatic feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getAromaticFeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, const Vector3D& featLoc,
    double scale) {
    std::string dirType = "linear";
    Vector3D head = featLoc;

    if (static_cast<int>(featAtoms.size()) < 2) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), dirType);
    }

    // Get positions of feature atoms
    auto pos0 = conf.getAtomPos(featAtoms[0]);
    auto pos1 = conf.getAtomPos(featAtoms[1]);
    Vector3D p0(pos0.x, pos0.y, pos0.z);
    Vector3D p1(pos1.x, pos1.y, pos1.z);

    Vector3D v1 = p0 - head;
    Vector3D v2 = p1 - head;
    Vector3D norm1 = v1.cross(v2);
    norm1.normalize();
    norm1 = norm1 * scale;

    Vector3D norm2 = head - norm1;
    norm1 = norm1 + head;

    std::vector<std::pair<Vector3D, Vector3D>> result;
    result.push_back(std::make_pair(head, norm1));
    result.push_back(std::make_pair(head, norm2));

    return std::make_pair(result, dirType);
}

// Get acceptor 2 feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getAcceptor2FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    int aid = featAtoms[0];
    auto cpt_pos = conf.getAtomPos(aid);
    Vector3D cpt(cpt_pos.x, cpt_pos.y, cpt_pos.z);

    const auto* mol = &conf.getOwningMol();
    // auto* atom = mol->getAtomWithIdx(aid); // Unused variable

    std::vector<const RDKit::Atom*> hydrogens;
    std::vector<const RDKit::Atom*> heavy_nbrs;

    // Find neighbors by checking bonds - RDKit C++ API doesn't have getNeighbors()
    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
        if (static_cast<int>(j) != aid) {
            auto bond = mol->getBondBetweenAtoms(aid, j);
            if (bond) {
                auto* nbr = mol->getAtomWithIdx(j);
                if (nbr->getAtomicNum() == 1) {
                    hydrogens.push_back(nbr);
                } else {
                    heavy_nbrs.push_back(nbr);
                }
            }
        }
    }

    if (heavy_nbrs.size() != 2) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    Vector3D bvec = findAvgVec(conf, cpt, heavy_nbrs);
    bvec = bvec * (-1.0 * scale);

    if (mol->getAtomWithIdx(aid)->getAtomicNum() == 8) {  // Oxygen
        // Assume sp3 hybridization
        auto pos0 = conf.getAtomPos(heavy_nbrs[0]->getIdx());
        auto pos1 = conf.getAtomPos(heavy_nbrs[1]->getIdx());
        Vector3D v1(pos0.x, pos0.y, pos0.z);
        Vector3D v2(pos1.x, pos1.y, pos1.z);

        v1 = v1 - cpt;
        v2 = v2 - cpt;
        Vector3D rotAxis = v1 - v2;
        rotAxis.normalize();

        Vector3D bv1 = arbAxisRotation(54.5, rotAxis, bvec);
        bv1 = bv1 + cpt;

        Vector3D bv2 = arbAxisRotation(-54.5, rotAxis, bvec);
        bv2 = bv2 + cpt;

        std::vector<std::pair<Vector3D, Vector3D>> result;
        result.push_back(std::make_pair(cpt, bv1));
        result.push_back(std::make_pair(cpt, bv2));

        return std::make_pair(result, "linear");
    } else {
        bvec = bvec + cpt;
        std::vector<std::pair<Vector3D, Vector3D>> result;
        result.push_back(std::make_pair(cpt, bvec));
        return std::make_pair(result, "linear");
    }
}

// Check planarity
bool checkPlanarity(const RDKit::Conformer& conf, const Vector3D& cpt,
                    const std::vector<const RDKit::Atom*>& nbrs, double tol) {
    if (nbrs.size() != 3) {
        return false;
    }

    auto pos0 = conf.getAtomPos(nbrs[0]->getIdx());
    auto pos1 = conf.getAtomPos(nbrs[1]->getIdx());
    auto pos2 = conf.getAtomPos(nbrs[2]->getIdx());

    Vector3D v1(pos0.x, pos0.y, pos0.z);
    Vector3D v2(pos1.x, pos1.y, pos1.z);
    Vector3D v3(pos2.x, pos2.y, pos2.z);

    v1 = v1 - cpt;
    v2 = v2 - cpt;
    v3 = v3 - cpt;

    Vector3D normal = v1.cross(v2);
    double dotP = std::abs(v3.dot(normal));

    return dotP <= tol;
}

// Get tetrahedral feature vector
std::vector<std::pair<Vector3D, Vector3D>> getTetrahedralFeatVect(const RDKit::Conformer& conf,
                                                                  int aid, double scale) {
    const auto* mol = &conf.getOwningMol();
    auto cpt_pos = conf.getAtomPos(aid);
    Vector3D cpt(cpt_pos.x, cpt_pos.y, cpt_pos.z);

    // auto* atom = mol->getAtomWithIdx(aid); // Unused variable
    // Find neighbors by checking bonds - RDKit C++ API doesn't have getNeighbors()
    std::vector<const RDKit::Atom*> nbrs;
    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
        if (static_cast<int>(j) != aid) {
            auto bond = mol->getBondBetweenAtoms(aid, j);
            if (bond) {
                nbrs.push_back(mol->getAtomWithIdx(j));
            }
        }
    }

    std::vector<std::pair<Vector3D, Vector3D>> result;

    if (!checkPlanarity(conf, cpt, nbrs, 0.1)) {
        Vector3D bvec = findAvgVec(conf, cpt, nbrs);
        bvec = bvec * (-1.0 * scale);
        bvec = bvec + cpt;
        result.push_back(std::make_pair(cpt, bvec));
    }

    return result;
}

// Get donor 3 feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getDonor3FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    int aid = featAtoms[0];
    auto tfv = getTetrahedralFeatVect(conf, aid, scale);
    return std::make_pair(tfv, "linear");
}

// Get acceptor 3 feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getAcceptor3FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    int aid = featAtoms[0];
    auto tfv = getTetrahedralFeatVect(conf, aid, scale);
    return std::make_pair(tfv, "linear");
}

// Get donor 2 feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getDonor2FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    int aid = featAtoms[0];
    const auto* mol = &conf.getOwningMol();
    auto cpt_pos = conf.getAtomPos(aid);
    Vector3D cpt(cpt_pos.x, cpt_pos.y, cpt_pos.z);

    // auto* atom = mol->getAtomWithIdx(aid); // Unused variable
    // Find neighbors by checking bonds - RDKit C++ API doesn't have getNeighbors()
    std::vector<const RDKit::Atom*> nbrs;
    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
        if (static_cast<int>(j) != aid) {
            auto bond = mol->getBondBetweenAtoms(aid, j);
            if (bond) {
                nbrs.push_back(mol->getAtomWithIdx(j));
            }
        }
    }

    if (nbrs.size() < 2) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "linear");
    }

    std::vector<const RDKit::Atom*> hydrogens;
    std::vector<const RDKit::Atom*> heavy_nbrs;

    for (auto* nbr : nbrs) {
        if (nbr->getAtomicNum() == 1) {
            hydrogens.push_back(nbr);
        } else {
            heavy_nbrs.push_back(nbr);
        }
    }

    std::vector<std::pair<Vector3D, Vector3D>> result;

    if (heavy_nbrs.size() == 2) {
        // No hydrogens expected
        Vector3D bvec = findAvgVec(conf, cpt, heavy_nbrs);
        bvec = bvec * (-1.0 * scale);
        bvec = bvec + cpt;
        result.push_back(std::make_pair(cpt, bvec));
    } else if (heavy_nbrs.size() == 3 && hydrogens.size() == 1) {
        // Check if hydrogen is in plane (sp2) or out of plane (sp3)
        int hid = hydrogens[0]->getIdx();
        auto h_pos = conf.getAtomPos(hid);
        Vector3D bvec(h_pos.x, h_pos.y, h_pos.z);
        bvec = bvec - cpt;
        bvec.normalize();
        bvec = bvec * scale;
        bvec = bvec + cpt;

        if (checkPlanarity(conf, cpt, heavy_nbrs)) {
            // Only hydrogen direction needed
            result.push_back(std::make_pair(cpt, bvec));
        } else {
            // Non-planar configuration - assume sp3
            Vector3D ovec = findAvgVec(conf, cpt, heavy_nbrs);
            ovec = ovec * (-1.0 * scale);
            ovec = ovec + cpt;
            result.push_back(std::make_pair(cpt, bvec));
            result.push_back(std::make_pair(cpt, ovec));
        }
    } else if (heavy_nbrs.size() >= 4) {
        // Two or more hydrogens - use their directions
        for (auto* hydrogen : hydrogens) {
            int hid = hydrogen->getIdx();
            auto h_pos = conf.getAtomPos(hid);
            Vector3D bvec(h_pos.x, h_pos.y, h_pos.z);
            bvec = bvec - cpt;
            bvec.normalize();
            bvec = bvec * scale;
            bvec = bvec + cpt;
            result.push_back(std::make_pair(cpt, bvec));
        }
    }

    return std::make_pair(result, "linear");
}

// Get donor 1 feature vectors
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getDonor1FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "cone");
    }

    int aid = featAtoms[0];
    const auto* mol = &conf.getOwningMol();
    // auto* atom = mol->getAtomWithIdx(aid); // Unused variable
    // Find neighbors by checking bonds - RDKit C++ API doesn't have getNeighbors()
    std::vector<const RDKit::Atom*> nbrs;
    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
        if (static_cast<int>(j) != aid) {
            auto bond = mol->getBondBetweenAtoms(aid, j);
            if (bond) {
                nbrs.push_back(mol->getAtomWithIdx(j));
            }
        }
    }

    // Find the neighboring heavy atom
    int hnbr = -1;
    for (auto* nbr : nbrs) {
        if (nbr->getAtomicNum() != 1) {
            hnbr = nbr->getIdx();
            break;
        }
    }

    if (hnbr == -1) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "cone");
    }

    auto cpt_pos = conf.getAtomPos(aid);
    auto hnbr_pos = conf.getAtomPos(hnbr);
    Vector3D cpt(cpt_pos.x, cpt_pos.y, cpt_pos.z);
    Vector3D v1(hnbr_pos.x, hnbr_pos.y, hnbr_pos.z);

    v1 = v1 - cpt;
    v1.normalize();
    v1 = v1 * (-1.0 * scale);
    v1 = v1 + cpt;

    std::vector<std::pair<Vector3D, Vector3D>> result;
    result.push_back(std::make_pair(cpt, v1));

    return std::make_pair(result, "cone");
}

// Get acceptor 1 feature vectors - this is the most complex one
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> getAcceptor1FeatVects(
    const RDKit::Conformer& conf, const std::vector<int>& featAtoms, double scale) {
    if (featAtoms.empty()) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "cone");
    }

    int aid = featAtoms[0];
    const auto* mol = &conf.getOwningMol();
    // auto* atom = mol->getAtomWithIdx(aid); // Unused variable
    // Find neighbors by checking bonds - RDKit C++ API doesn't have getNeighbors()
    std::vector<const RDKit::Atom*> nbrs;
    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
        if (static_cast<int>(j) != aid) {
            auto bond = mol->getBondBetweenAtoms(aid, j);
            if (bond) {
                nbrs.push_back(mol->getAtomWithIdx(j));
            }
        }
    }

    auto cpt_pos = conf.getAtomPos(aid);
    Vector3D cpt(cpt_pos.x, cpt_pos.y, cpt_pos.z);

    // Find the adjacent heavy atom
    const RDKit::Atom* heavyAt = nullptr;
    for (auto* nbr : nbrs) {
        if (nbr->getAtomicNum() != 1) {
            heavyAt = nbr;
            break;
        }
    }

    if (heavyAt == nullptr) {
        return std::make_pair(std::vector<std::pair<Vector3D, Vector3D>>(), "cone");
    }

    // Check bond type
    auto bond = mol->getBondBetweenAtoms(aid, heavyAt->getIdx());
    bool singleBnd = bond && bond->getBondType() == RDKit::Bond::SINGLE;

    // Special case for sulfur
    bool sulfur = heavyAt->getAtomicNum() == 16;

    if (singleBnd || sulfur) {
        // Use bond direction (inverted)
        auto heavy_pos = conf.getAtomPos(heavyAt->getIdx());
        Vector3D v1(heavy_pos.x, heavy_pos.y, heavy_pos.z);
        v1 = v1 - cpt;
        v1.normalize();
        v1 = v1 * (-1.0 * scale);
        v1 = v1 + cpt;

        std::vector<std::pair<Vector3D, Vector3D>> result;
        result.push_back(std::make_pair(cpt, v1));
        return std::make_pair(result, "cone");
    } else {
        // Assume sp2 hybridization and create two direction vectors
        auto heavy_pos = conf.getAtomPos(heavyAt->getIdx());
        Vector3D v1(heavy_pos.x, heavy_pos.y, heavy_pos.z);

        // Find another neighbor for plane definition
        Vector3D plane_vector(1.0, 0.0, 0.0);  // Default
        // Find heavy atom neighbors by checking bonds
        std::vector<const RDKit::Atom*> heavy_nbrs;
        for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
            if (j != heavyAt->getIdx()) {
                auto bond = mol->getBondBetweenAtoms(heavyAt->getIdx(), j);
                if (bond) {
                    heavy_nbrs.push_back(mol->getAtomWithIdx(j));
                }
            }
        }
        for (auto* nbr : heavy_nbrs) {
            if (static_cast<int>(nbr->getIdx()) != aid) {
                auto nbr_pos = conf.getAtomPos(nbr->getIdx());
                plane_vector = Vector3D(nbr_pos.x, nbr_pos.y, nbr_pos.z) - v1;
                break;
            }
        }

        Vector3D pt1 = plane_vector;
        pt1 = pt1 - v1;
        v1 = v1 - cpt;

        Vector3D rotAxis = v1.cross(pt1);
        rotAxis.normalize();

        Vector3D bv1 = arbAxisRotation(120, rotAxis, v1);
        bv1.normalize();
        bv1 = bv1 * scale;
        bv1 = bv1 + cpt;

        Vector3D bv2 = arbAxisRotation(-120, rotAxis, v1);
        bv2.normalize();
        bv2 = bv2 * scale;
        bv2 = bv2 + cpt;

        std::vector<std::pair<Vector3D, Vector3D>> result;
        result.push_back(std::make_pair(cpt, bv1));
        result.push_back(std::make_pair(cpt, bv2));

        return std::make_pair(result, "linear");
    }
}

// Find hydrogen atoms
std::vector<int> findHydAtoms(const std::vector<int>& nbrs,
                              const std::vector<std::string>& atomNames) {
    std::vector<int> hAtoms;
    for (int nid : nbrs) {
        if (nid >= 0 && nid < static_cast<int>(atomNames.size()) && atomNames[nid] == "H") {
            hAtoms.push_back(nid);
        }
    }
    return hAtoms;
}

}  // namespace RDKitUtils
