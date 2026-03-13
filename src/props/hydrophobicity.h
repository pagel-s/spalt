/**
 * @file hydrophobicity.h
 * @brief Hydrophobicity surface property computation
 * @author Sebastian
 * @date 2025
 */

#ifndef HYDROPHOBICITY_H
#define HYDROPHOBICITY_H

#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/GraphMol.h>
#include <map>
#include <string>
#include <vector>
#include "surface_property.h"

/**
 * @class HydrophobicityProperty
 * @brief Computes hydrophobicity values for surface vertices using fragment-based LogP calculations
 *
 * This class implements hydrophobicity computation based on fragment-based LogP calculations.
 * It uses RDKit's BRICS fragmentation and Crippen LogP to assign hydrophobicity values
 * to each atom based on the fragment it belongs to.
 *
 * The hydrophobicity values are computed using the following approach:
 * 1. Fragment the molecule using BRICS bonds
 * 2. Calculate LogP for each fragment using Crippen method
 * 3. Convert LogP to hydrophobicity using kd_from_logp function
 * 4. Assign hydrophobicity values to surface vertices based on atom names
 */
class HydrophobicityProperty : public ISurfaceProperty {
  public:
    /**
     * @brief Construct a new HydrophobicityProperty object
     */
    HydrophobicityProperty();

    /**
     * @brief Destroy the HydrophobicityProperty object
     */
    virtual ~HydrophobicityProperty();

    /**
     * @brief Compute hydrophobicity values for surface vertices
     *
     * @param surface The surface object containing vertices
     * @param cache Property cache for storing computed values
     */
    virtual void compute(const Surface& surface,
                         std::unordered_map<std::string, std::any>& cache) override;

    /**
     * @brief Get the property key
     *
     * @return Property key string
     */
    virtual std::string key() const override {
        return "hydrophobicity";
    }

    // Public methods for testing
    /**
     * @brief Convert LogP value to hydrophobicity (kd) value
     *
     * @param logp LogP value from Crippen calculation
     * @param kd_min Minimum kd value (default: -4.5)
     * @param kd_max Maximum kd value (default: 4.5)
     * @return Hydrophobicity (kd) value
     */
    double kdFromLogP(double logp, double kd_min = -4.5, double kd_max = 4.5) const;

    /**
     * @brief Assign hydrophobicity values to atoms based on fragment analysis
     *
     * @param mol RDKit molecule object
     * @return Map from atom names to hydrophobicity values
     */
    std::map<std::string, double> molAssignKd(std::shared_ptr<RDKit::ROMol> mol) const;

  private:
    /**
     * @brief Get molecular fragments using BRICS fragmentation
     *
     * @param mol RDKit molecule object
     * @return Vector of molecular fragments
     */
    std::vector<RDKit::ROMOL_SPTR> getFragments(std::shared_ptr<RDKit::ROMol> mol) const;
};

#endif  // HYDROPHOBICITY_H
