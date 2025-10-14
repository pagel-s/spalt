/**
 * @file esp.h
 * @brief Electrostatic Potential (ESP) property computation for molecular surfaces
 * @author Sebastian
 * @date 2024
 */

#pragma once
#include <vector>
#include <string>
#include "surface_property.h"

/**
 * @enum ChargeMethod
 * @brief Enumeration of available charge calculation methods
 * 
 * Defines the different methods available for computing atomic charges
 * used in electrostatic potential calculations.
 */
enum class ChargeMethod {
    XTB,    ///< Use xtb program for charge calculation (external quantum chemistry)
    RDKIT   ///< Use RDKit built-in charge calculation methods (Gasteiger charges)
};

/**
 * @struct SurfaceESPProperty
 * @brief Surface property for computing Electrostatic Potential (ESP) charges
 * 
 * The SurfaceESPProperty class computes atomic charges for electrostatic potential
 * calculations using either external quantum chemistry programs (xtb) or built-in
 * RDKit methods (Gasteiger charges). It supports both methods through a unified
 * interface with automatic fallback mechanisms.
 * 
 * Key features:
 * - Dual charge calculation methods (XTB and RDKit)
 * - Constructor-based method configuration
 * - Temporary directory management for external programs
 * - Fail-safe cleanup of temporary files
 * - Caching of computed results
 * 
 * @see ISurfaceProperty
 * @see Surface
 */
struct SurfaceESPProperty : ISurfaceProperty {
  /**
   * @brief Construct a new SurfaceESPProperty object
   * 
   * @param method Charge calculation method to use (default: ChargeMethod::XTB)
   */
  SurfaceESPProperty(ChargeMethod method = ChargeMethod::XTB) : charge_method_(method) {}
  
  /**
   * @brief Get the property key identifier
   * 
   * @return String identifier "surface.esp"
   */
  std::string key() const override { return "surface.esp"; }
  
  /**
   * @brief Compute ESP charges using the configured method
   * 
   * Computes atomic charges using either XTB or RDKit based on the method
   * configured in the constructor. Results are stored in the cache.
   * 
   * @param s Surface object containing the molecular structure
   * @param cache Property cache for storing computed results
   * 
   * @throws std::runtime_error if charge computation fails
   */
  void compute(const Surface& s, std::unordered_map<std::string, std::any>& cache) override;
  
  /**
   * @brief Legacy compute method with explicit charge method parameter
   * 
   * @deprecated Use constructor parameter instead: SurfaceESPProperty(ChargeMethod::RDKIT)
   * 
   * @param s Surface object containing the molecular structure
   * @param cache Property cache for storing computed results
   * @param method Charge calculation method to use
   * 
   * @throws std::runtime_error if charge computation fails
   */
  [[deprecated("Use constructor parameter instead: SurfaceESPProperty(ChargeMethod::RDKIT)")]]
  void compute(const Surface& s, std::unordered_map<std::string, std::any>& cache, ChargeMethod method);
  
  /**
   * @brief Load ESP charges from file
   * 
   * Reads atomic charges from a text file where each line contains a single charge value.
   * Used primarily for loading XTB-computed charges from output files.
   * 
   * @param esp_file Path to the ESP charges file (default: "xtb_esp.dat")
   * @return Vector of atomic charges
   * 
   * @throws std::runtime_error if file cannot be opened or contains invalid data
   */
  std::vector<double> load_esp(std::string_view esp_file = "xtb_esp.dat") const;
      
  /**
   * @brief Static helper to retrieve ESP charges from property cache
   * 
   * @param cache Property cache containing computed ESP charges
   * @return Vector of atomic charges
   * 
   * @throws std::runtime_error if ESP charges not found in cache
   * @throws std::bad_any_cast if cache contains wrong data type
   */
  static std::vector<double> get_esp_charges(const std::unordered_map<std::string, std::any>& cache);
  
  /**
   * @brief Compute electrostatic potential values at surface points
   * 
   * Calculates the Coulomb potential at each surface vertex using the computed atomic charges
   * and atomic positions. This implements the electrostatics calculation given point charges.
   * 
   * @param charges Vector of atomic charges (from compute_xtb_charges or compute_rdkit_charges)
   * @param positions Vector of atomic positions (x,y,z coordinates for each atom)
   * @param surface_points Vector of surface point coordinates (x,y,z coordinates for each vertex)
   * @return Vector of electrostatic potential values at each surface point
   * 
   * @throws std::runtime_error if input vectors have mismatched dimensions
   */
  static std::vector<double> compute_electrostatic_potential(
      const std::vector<double>& charges,
      const std::vector<std::vector<double>>& positions,
      const std::vector<std::vector<double>>& surface_points);
  
  /**
   * @brief Get atomic positions from a surface's molecule
   * 
   * Extracts the 3D coordinates of all atoms from the molecule associated with the surface.
   * 
   * @param surface Surface object containing the molecule
   * @return Vector of atomic positions (x,y,z coordinates for each atom)
   * 
   * @throws std::runtime_error if molecule has no conformers or conformer access fails
   */
  static std::vector<std::vector<double>> get_atomic_positions(const Surface& surface);
  
  /**
   * @brief Get surface point coordinates from a surface
   * 
   * Extracts the 3D coordinates of all surface vertices.
   * 
   * @param surface Surface object
   * @return Vector of surface point coordinates (x,y,z coordinates for each vertex)
   */
  static std::vector<std::vector<double>> get_surface_point_coordinates(const Surface& surface);
  
  /**
   * @brief Compute ESP charges and calculate electrostatic potentials at surface points
   * 
   * This is a convenience method that combines charge computation and ESP calculation.
   * It computes atomic charges using the configured method, then calculates the electrostatic
   * potential at each surface vertex.
   * 
   * @param surface Surface object containing the molecular structure
   * @param cache Property cache for storing computed results
   * @return Vector of electrostatic potential values at each surface point
   * 
   * @throws std::runtime_error if charge computation or ESP calculation fails
   */
  std::vector<double> compute_esp_potentials(const Surface& surface, std::unordered_map<std::string, std::any>& cache);
  
  /**
   * @brief Static helper to retrieve ESP potentials from property cache
   * 
   * @param cache Property cache containing computed ESP potentials
   * @return Vector of electrostatic potential values at surface points
   * 
   * @throws std::runtime_error if ESP potentials not found in cache
   * @throws std::bad_any_cast if cache contains wrong data type
   */
  static std::vector<double> get_esp_potentials(const std::unordered_map<std::string, std::any>& cache);
  
  /**
   * @brief Get the configured charge calculation method
   * 
   * @return The charge method configured in the constructor
   */
  ChargeMethod getChargeMethod() const { return charge_method_; }
  
private:
  /**
   * @brief Configured charge calculation method
   */
  ChargeMethod charge_method_;
  
  /**
   * @brief Compute charges using XTB program
   * 
   * Runs the xtb quantum chemistry program in a temporary directory to compute
   * atomic charges. Handles temporary directory creation, cleanup, and error recovery.
   * 
   * @param s Surface object containing the molecular structure
   * @return Vector of computed atomic charges
   * 
   * @throws std::runtime_error if XTB computation fails or temporary directory issues
   */
  std::vector<double> compute_xtb_charges(const Surface& s);
  
  /**
   * @brief Compute charges using RDKit built-in methods
   * 
   * Uses RDKit's Gasteiger charge calculation with fallback to formal charges.
   * This method is faster and doesn't require external programs but may be less accurate.
   * 
   * @param s Surface object containing the molecular structure
   * @param cache Property cache for storing computation metadata
   * @return Vector of computed atomic charges
   * 
   * @throws std::runtime_error if RDKit charge computation fails
   */
  std::vector<double> compute_rdkit_charges(const Surface& s, std::unordered_map<std::string, std::any>& cache);
};