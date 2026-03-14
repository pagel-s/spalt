/**
 * @file pharmacophore.h
 * @brief Pharmacophore surface property computation
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>
#include "surface_property.h"

namespace RDKit {
class ROMol;
}

enum class PharmaType { ACCEPTOR, AROMATIC, DONOR, HYDROPHOBE, LUMPED_HYDROPHOBE, NEG_IONIZABLE, POS_IONIZABLE, ZN_BINDER };

class PharmacophoreProperty : public ISurfaceProperty {
  public:
    PharmacophoreProperty(PharmaType type);
    virtual ~PharmacophoreProperty() = default;

    std::string key() const override;
    std::string getFamilyName() const;
    void compute(const Surface& surface, std::unordered_map<std::string, std::any>& cache) override;

  private:
    PharmaType type_;

    struct Feature {
        Eigen::Vector3d position;
        bool has_direction;
        Eigen::Vector3d direction;  // Normalized direction vector
    };

    std::vector<Feature> extractFeatures(std::shared_ptr<RDKit::ROMol> mol) const;
};
