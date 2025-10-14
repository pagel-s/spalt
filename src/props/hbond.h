/**
 * @file hbond.h
 * @brief Hydrogen bond potential surface property computation
 * @author Sebastian
 * @date 2025
 */

#ifndef HBOND_H
#define HBOND_H

#include "surface_property.h"
#include "rdkit_utils.h"
#include <map>
#include <vector>
#include <string>
#include <set>
#include <Eigen/Core>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureDef.h>

/**
 * @class HBondProperty
 * @brief Computes hydrogen bond potential values for surface vertices
 * 
 * This class implements hydrogen bond potential computation based on donor/acceptor analysis.
 * It uses RDKit's chemical features to identify donor and acceptor atoms, then calculates
 * hydrogen bond potentials based on geometric considerations and angle deviations.
 * 
 * The hydrogen bond potential computation follows this approach:
 * 1. Identify donor hydrogens and acceptor atoms using RDKit chemical features
 * 2. For donor hydrogens: calculate angle deviation from ideal 180° (π) direction
 * 3. For acceptor atoms: calculate angle deviation from ideal acceptor directions
 * 4. Apply angle penalty functions to get hydrogen bond potential values
 * 5. Assign values to surface vertices based on atom names
 */
class HBondProperty : public ISurfaceProperty {
public:
    /**
     * @brief Construct a new HBondProperty object
     */
    HBondProperty();
    
    /**
     * @brief Destroy the HBondProperty object
     */
    virtual ~HBondProperty();
    
    /**
     * @brief Compute hydrogen bond potential values for surface vertices
     * 
     * @param surface The surface object containing vertices
     * @param cache Property cache for storing computed values
     */
    virtual void compute(const Surface& surface, std::unordered_map<std::string, std::any>& cache) override;
    
    /**
     * @brief Get the property key
     * 
     * @return Property key string
     */
    virtual std::string key() const override { return "hbond"; }
    
    /**
     * @brief Get donor hydrogens from the molecule
     * 
     * @param mol RDKit molecule object
     * @return Set of donor hydrogen atom indices
     */
    std::set<int> getDonorHydrogens(std::shared_ptr<RDKit::ROMol> mol) const;
    
    /**
     * @brief Get acceptor atoms from the molecule
     * 
     * @param mol RDKit molecule object
     * @return Set of acceptor atom indices
     */
    std::set<int> getAcceptors(std::shared_ptr<RDKit::ROMol> mol) const;
    
    /**
     * @brief Get donors from the molecule
     * 
     * @param mol RDKit molecule object
     * @return Set of donor atom indices
     */
    std::set<int> getDonors(std::shared_ptr<RDKit::ROMol> mol) const;
    
    /**
     * @brief Get acceptor type (number of heavy atom neighbors)
     * 
     * @param acceptor_atom RDKit atom object
     * @return Number of heavy atom neighbors (1, 2, or 3)
     */
    int getAcceptorType(const RDKit::Atom& acceptor_atom) const;

private:
    
    /**
     * @brief Compute angle deviation from ideal value
     * 
     * @param a First point
     * @param b Second point (vertex)
     * @param c Third point
     * @param ideal_angle Ideal angle in radians
     * @return Absolute deviation from ideal angle
     */
    double computeAngleDeviation(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
                                const Eigen::Vector3d& c, double ideal_angle) const;
    
    /**
     * @brief Compute angle penalty based on deviation
     * 
     * @param angle_deviation Angle deviation in radians
     * @return Penalty value between 0 and 1
     */
    double computeAnglePenalty(double angle_deviation) const;
    
    /**
     * @brief Compute hydrogen bond potential for a specific atom and surface point
     * This is the core helper function that matches the Python computeChargeHelperMol
     * 
     * @param mol RDKit molecule object
     * @param atom_idx Atom index
     * @param donor_hydrogens Set of donor hydrogen indices
     * @param acceptors Set of acceptor atom indices
     * @param surface_point Surface point coordinates
     * @return Hydrogen bond potential value
     */
    double computeChargeHelperMol(std::shared_ptr<RDKit::ROMol> mol, unsigned int atom_idx,
                                const std::set<int>& donor_hydrogens, 
                                const std::set<int>& acceptors,
                                const Eigen::Vector3d& surface_point) const;
    
    /**
     * @brief Legacy function - now calls computeChargeHelperMol
     */
    double computeHBondPotential(std::shared_ptr<RDKit::ROMol> mol, int atom_idx,
                                const std::set<int>& donor_hydrogens, 
                                const std::set<int>& acceptors,
                                const Eigen::Vector3d& surface_point) const;
    
    /**
     * @brief Get acceptor direction vectors for type 1 acceptors
     * 
     * @param mol RDKit molecule object
     * @param atom_idx Acceptor atom index
     * @return Vector of direction vectors and type
     */
    std::pair<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>, std::string> 
    getAcceptorFeatVects(std::shared_ptr<RDKit::ROMol> mol, int atom_idx) const;
    
    /**
     * @brief Standard deviation for hydrogen bond angles
     */
    static constexpr double HBOND_STD_DEV = M_PI / 3.0; // 60 degrees
};

#endif // HBOND_H
