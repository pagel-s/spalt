/**
 * @file pharmacophore.h
 * @brief Pharmacophore surface property computation
 */

#pragma once

#include "surface_property.h"
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Core>

namespace RDKit {
    class ROMol;
}

enum class PharmaType {
    AROMATIC,
    POSITIVE,
    NEGATIVE,
    DONOR,
    ACCEPTOR
};

class PharmacophoreProperty : public ISurfaceProperty {
public:
    PharmacophoreProperty(PharmaType type);
    virtual ~PharmacophoreProperty() = default;

    std::string key() const override;
    void compute(const Surface& surface, std::unordered_map<std::string, std::any>& cache) override;

private:
    PharmaType type_;

    struct Feature {
        Eigen::Vector3d position;
        bool has_direction;
        Eigen::Vector3d direction; // Normalized direction vector
    };

    std::vector<Feature> extractFeatures(std::shared_ptr<RDKit::ROMol> mol) const;
};
