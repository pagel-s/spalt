/**
 * @file property_params.cpp
 * @brief Property parameter definitions implementation
 * @author Sebastian
 * @date 2025
 */

#include "property_params.h"

// Property parameter definitions with standardized min/max values
const std::unordered_map<std::string, PropertyParams> PropertyParamsRegistry::registry_ = {
    {"esp", PropertyParams("esp", -10.0, 10.0)},  // Electrostatic potential in eV
    {"hydrophobicity", PropertyParams("hydrophobicity", -4.5, 4.5)},  // Hydrophobicity scale
    {"hbond", PropertyParams("hbond", -1.0, 1.0)},                    // Hydrogen bond potential
    {"pharma_aromatic", PropertyParams("pharma_aromatic", 0.0, 1.0)},
    {"pharma_pos", PropertyParams("pharma_pos", 0.0, 1.0)},
    {"pharma_neg", PropertyParams("pharma_neg", 0.0, 1.0)},
    {"pharma_donor", PropertyParams("pharma_donor", 0.0, 1.0)},
    {"pharma_acceptor", PropertyParams("pharma_acceptor", 0.0, 1.0)},
};
