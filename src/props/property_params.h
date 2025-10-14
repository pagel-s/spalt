/**
 * @file property_params.h
 * @brief Property parameter definitions for consistent color mapping
 * @author Sebastian
 * @date 2025
 */

#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

/**
 * @struct PropertyParams
 * @brief Parameters for surface property computation and visualization
 * 
 * Contains min/max values for consistent color mapping across different molecules
 * and property-specific configuration parameters. These parameters ensure that
 * surface properties are visualized with consistent color scales across different
 * molecular structures for meaningful comparison.
 */
struct PropertyParams {
    double min_value = 0.0;  ///< Minimum value for color mapping and normalization (default: 0.0)
    double max_value = 1.0;  ///< Maximum value for color mapping and normalization (default: 1.0)
    std::string name;        ///< Property name identifier (e.g., "esp", "hydrophobicity", "hbond")
    
    PropertyParams(const std::string& prop_name, double min_val, double max_val) 
        : min_value(min_val), max_value(max_val), name(prop_name) {}
};

/**
 * @class PropertyParamsRegistry
 * @brief Registry for property parameter definitions
 * 
 * Provides a centralized way to manage property parameters including
 * min/max values for consistent color mapping across molecules.
 */
class PropertyParamsRegistry {
public:
    /**
     * @brief Get property parameters by name
     * 
     * @param property_name Name of the property (esp, hydrophobicity, hbond, etc.)
     * @return PropertyParams struct with min/max values
     * @throws std::runtime_error if property not found
     */
    static PropertyParams getParams(const std::string& property_name) {
        auto it = registry_.find(property_name);
        if (it == registry_.end()) {
            throw std::runtime_error("Unknown property: " + property_name + 
                                   ". Available properties: esp, hydrophobicity, hbond");
        }
        return it->second;
    }
    
    /**
     * @brief Check if a property is registered
     * 
     * @param property_name Name of the property
     * @return true if property is registered, false otherwise
     */
    static bool hasProperty(const std::string& property_name) {
        return registry_.find(property_name) != registry_.end();
    }
    
    /**
     * @brief Get all registered property names
     * 
     * @return Vector of property names
     */
    static std::vector<std::string> getPropertyNames() {
        std::vector<std::string> names;
        for (const auto& pair : registry_) {
            names.push_back(pair.first);
        }
        return names;
    }

private:
    /**
     * @brief Property parameter registry with predefined min/max values
     * 
     * Defines standard ranges for consistent color mapping:
     * - ESP: -10 to +10 eV (electrostatic potential)
     * - Hydrophobicity: -4.5 to +4.5 (hydrophobicity scale)
     * - HBond: -1 to +1 (hydrogen bond potential)
     */
    static const std::unordered_map<std::string, PropertyParams> registry_;
};

// Property parameter definitions are in property_params.cpp
