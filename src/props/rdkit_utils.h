/**
 * @file rdkit_utils.h
 * @brief RDKit utility functions for chemical feature calculations
 * @author Sebastian
 * @date 2025
 * 
 * This file contains C++ implementations of RDKit's chemical feature utility functions,
 * converted from the Python feature factory code. These functions are used for
 * hydrogen bond potential calculations and other chemical feature computations.
 */

#ifndef RDKIT_UTILS_H
#define RDKIT_UTILS_H

#include <vector>
#include <string>
#include <memory>
#include <Eigen/Core>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Conformer.h>

namespace RDKitUtils {

/**
 * @brief 3D vector operations for chemical feature calculations
 */
class Vector3D {
public:
    double x, y, z;
    
    Vector3D() : x(0.0), y(0.0), z(0.0) {}
    Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vector3D(const Eigen::Vector3d& vec) : x(vec[0]), y(vec[1]), z(vec[2]) {}
    
    // Convert to Eigen vector
    Eigen::Vector3d toEigen() const { return Eigen::Vector3d(x, y, z); }
    
    // Vector operations
    Vector3D operator+(const Vector3D& other) const;
    Vector3D operator-(const Vector3D& other) const;
    Vector3D operator*(double scalar) const;
    Vector3D operator/(double scalar) const;
    Vector3D operator-() const;
    
    // Dot product
    double dot(const Vector3D& other) const;
    
    // Cross product
    Vector3D cross(const Vector3D& other) const;
    
    // Normalize vector
    void normalize();
    Vector3D normalized() const;
    
    // Vector length
    double norm() const;
    double normsq() const;
    
    // Angle between vectors
    double angle(const Vector3D& other) const;
};

/**
 * @brief Cross product of two vectors
 */
Vector3D cross(const Vector3D& v1, const Vector3D& v2);

/**
 * @brief Find neighboring atom IDs from adjacency matrix
 * @param atomId Atom of interest
 * @param adjMat Adjacency matrix for the compound
 * @return Vector of neighboring atom IDs
 */
std::vector<int> findNeighbors(int atomId, const std::vector<std::vector<int>>& adjMat);

/**
 * @brief Find the average of normalized vectors going from center atom to neighbors
 * @param conf Conformer
 * @param center Center atom position
 * @param nbrs Neighbor atoms
 * @return Average normalized vector
 */
Vector3D findAvgVec(const RDKit::Conformer& conf, const Vector3D& center, 
                    const std::vector<const RDKit::Atom*>& nbrs);

/**
 * @brief Arbitrary axis rotation
 * @param theta Rotation angle in degrees
 * @param ax Rotation axis
 * @param pt Point to rotate
 * @return Rotated point
 */
Vector3D arbAxisRotation(double theta, const Vector3D& ax, const Vector3D& pt);

/**
 * @brief Get direction vectors for aromatic feature
 * @param conf Conformer
 * @param featAtoms List of atom IDs that make up the feature
 * @param featLoc Location of the aromatic feature
 * @param scale Size of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getAromaticFeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                    const Vector3D& featLoc, double scale = 1.5);

/**
 * @brief Get direction vectors for Acceptor of type 2
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getAcceptor2FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                     double scale = 1.5);

/**
 * @brief Get direction vectors for Donor of type 3
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getDonor3FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                  double scale = 1.5);

/**
 * @brief Get direction vectors for Acceptor of type 3
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getAcceptor3FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                     double scale = 1.5);

/**
 * @brief Get direction vectors for Donor of type 2
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getDonor2FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                  double scale = 1.5);

/**
 * @brief Get direction vectors for Donor of type 1
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getDonor1FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                  double scale = 1.5);

/**
 * @brief Get direction vectors for Acceptor of type 1
 * @param conf Conformer
 * @param featAtoms List of atoms that are part of the feature
 * @param scale Length of the direction vector
 * @return Pair of direction vectors and type
 */
std::pair<std::vector<std::pair<Vector3D, Vector3D>>, std::string> 
getAcceptor1FeatVects(const RDKit::Conformer& conf, const std::vector<int>& featAtoms, 
                     double scale = 1.5);

/**
 * @brief Check if atoms are planar
 * @param conf Conformer
 * @param cpt Center point
 * @param nbrs Neighbor atoms
 * @param tol Tolerance for planarity check
 * @return True if atoms are planar
 */
bool checkPlanarity(const RDKit::Conformer& conf, const Vector3D& cpt, 
                   const std::vector<const RDKit::Atom*>& nbrs, double tol = 1.0e-3);

/**
 * @brief Get tetrahedral feature vector
 * @param conf Conformer
 * @param aid Atom ID
 * @param scale Scale factor
 * @return Tetrahedral feature vector
 */
std::vector<std::pair<Vector3D, Vector3D>> 
getTetrahedralFeatVect(const RDKit::Conformer& conf, int aid, double scale);

/**
 * @brief Find hydrogen atoms in neighbor list
 * @param nbrs Neighbor atom IDs
 * @param atomNames Atom names array
 * @return Vector of hydrogen atom IDs
 */
std::vector<int> findHydAtoms(const std::vector<int>& nbrs, 
                             const std::vector<std::string>& atomNames);

} // namespace RDKitUtils

#endif // RDKIT_UTILS_H
