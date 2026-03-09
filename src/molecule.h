/**
 * @file molecule.h
 * @brief Molecular structure representation and manipulation using RDKit
 * @author Sebastian
 * @date 2024
 */

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <functional>
#include <GraphMol/GraphMol.h>
#include <Eigen/Geometry>
#include "conformer_params.h"
#include "props/esp.h"

// Forward declarations
class Surface;
class ISurfaceProperty;
class ConformerGenerator;

/**
 * @struct SurfaceParams
 * @brief Parameters for surface generation and molecular property computation
 * 
 * Comprehensive parameter structure controlling all aspects of molecular surface
 * generation, property computation, and output formatting. Used throughout the
 * SPALT system for consistent surface generation across different algorithms.
 */
struct SurfaceParams {
    int num_vertices = 1000;  ///< Number of vertices to generate on the molecular surface (default: 1000)
    double radius = 1.0;      ///< Probe radius for surface generation in Angstroms (default: 1.2)
    double density = 3.0;     ///< Surface density parameter for MSMS algorithm (default: 3.0)
    double hdensity = 3.0;    ///< High density parameter for MSMS algorithm (default: 3.0)
    std::string type = "tses"; ///< Surface type: "tses" (triangulated solvent excluded surface) or "ases" (analytical surface) (default: "tses")
    std::string program = "msms"; ///< Surface generation program: "msms" or "fibonacci" (default: "msms")
    std::string sample_method = "full"; ///< Vertex sampling method: "full" (all vertices), "fps" (farthest point sampling), or "random" (default: "full")
    std::vector<std::string> properties = {"esp", "hydrophobicity", "hbond"}; ///< List of surface properties to compute (default: compute all available properties)
    std::string output_type = "mesh"; ///< Output format: "mesh" (triangular mesh) or "points" (point cloud) (default: "mesh")
    bool include_normals = false; ///< Whether to include vertex normals in output (default: false for better PyMOL compatibility)
    int top_n = -1;  ///< Number of top conformers to save (-1 = save all, >0 = save only top N conformers) (default: -1)
    
    // Advanced conformer generation parameters
    ConformerParams conformer_params; ///< Advanced conformer generation settings including clustering and optimization parameters
    bool save_meshes = false; ///< Whether to save individual mesh files for each conformer (default: false)
    
    // ESP charge calculation method
    ChargeMethod charge_method = ChargeMethod::RDKIT; ///< Charge calculation method for ESP: XTB (quantum chemistry) or RDKIT (Gasteiger charges) (default: RDKIT)
    
    // Hydrogen control parameters
    bool addH = true;  ///< Whether to add explicit hydrogens to molecules (default: true)
};

/**
 * @class Molecule
 * @brief Wrapper class for RDKit molecular structure representation
 * 
 * The Molecule class provides a high-level interface for working with molecular
 * structures using RDKit. It handles SMILES parsing, conformer generation,
 * charge calculation, and file I/O operations.
 * 
 * Key features:
 * - SMILES string parsing and validation
 * - 3D conformer generation using RDKit
 * - Formal charge calculation
 * - Molecular file format conversion (molfile)
 * - Surface generation and property computation
 * - Extensible property system
 * 
 * @see Surface
 */
class Molecule {
    public:
        /**
         * @brief Construct a new Molecule object from SMILES string
         * 
         * Parses the SMILES string, validates the molecular structure, calculates
         * formal charge, and optionally generates conformers for 3D representation.
         * 
         * @param smiles SMILES string representation of the molecule
         * 
         * @throws std::runtime_error if SMILES parsing fails or molecule is invalid
         */
        Molecule(const char* smiles);
        
        /**
         * @brief Construct a new Molecule object from a molecular file
         * 
         * Supports SDF file format and uses existing 3D coordinates when available
         * 
         * @param filename Path to the molecular file
         * 
         * @throws std::runtime_error if file cannot be read or parsed
         */
        Molecule(const std::string& filename);
        
        /**
         * @brief Move constructor
         */
        Molecule(Molecule&& other) noexcept;
        
        /**
         * @brief Move assignment operator
         */
        Molecule& operator=(Molecule&& other) noexcept;
        
        /**
         * @brief Destroy the Molecule object
         */
        ~Molecule();

        /**
         * @brief Get the original SMILES string
         * 
         * @return The SMILES string used to construct the molecule
         */
        std::string get_smiles() const;
        
        /**
         * @brief Get a copy of the RDKit molecule object
         * 
         * @return Unique pointer to a copy of the RDKit ROMol object
         * @return nullptr if the molecule is invalid
         */
        std::unique_ptr<RDKit::ROMol> get_mol() const;

        /**
         * @brief Convert molecule to molfile format and save to file
         * 
         * Writes the molecular structure in MDL molfile format to the specified file.
         * This is commonly used for interfacing with other molecular modeling programs.
         * 
         * @param filename Output file path (default: "mol.mol")
         * @return true if successful, false otherwise
         * 
         * @note Returns false if molecule is invalid or file writing fails
         */
        bool mol_to_molfile(std::string_view filename = "mol.mol") const;

        /**
         * @brief Generate 3D conformers for the molecule
         * 
         * Uses RDKit's distance geometry algorithm to generate multiple 3D conformations
         * of the molecule. This is essential for surface generation and molecular visualization.
         * 
         * @param num_conformers Number of conformers to generate (default: 10)
         * 
         * @note Only generates conformers if the molecule is valid
         * @note May generate fewer conformers than requested if the algorithm fails
         */
        void calculate_conformers(unsigned int num_conformers = 10) const;
        
        /**
         * @brief Generate conformers using advanced clustering-based approach
         * 
         * Uses multiple sampling strategies (optimized, ETKDG, uniform) with K-means clustering
         * to select diverse, high-quality conformers. Based on the Python reference implementation.
         * 
         * @param params Conformer generation parameters (default: ConformerParams())
         * @return Map of cluster_id -> conformer_id for selected conformers
         * 
         * @note This method replaces existing conformers by default
         * @note Uses MMFF optimization, ETKDG sampling, and uniform random sampling
         * @note Applies K-means clustering to select diverse conformers
         */
        std::unordered_map<int, int> generateAdvancedConformers(const ConformerParams& params = ConformerParams());

        /**
         * @brief Get the formal charge of the molecule
         * 
         * @return Formal charge of the molecule
         */
        int get_charge() const;
        
        /**
         * @brief Compute a surface property for a specific conformer
         * 
         * Extensible API that automatically creates surface, registers property, and computes it.
         * Properties are registered through a factory system for easy extensibility.
         * 
         * @param property_name Name of the property to compute (e.g., "esp", "hydrophobicity")
         * @param conformer_id ID of the conformer (default: 0)
         * @param surface_params Surface generation parameters (default: SurfaceParams())
         * 
         * @throws std::runtime_error if property computation fails or property not registered
         * 
         * @example
         *   molecule.compute("esp");  // Compute ESP for conformer 0
         *   molecule.compute("esp", 1);  // Compute ESP for conformer 1
         */
        void compute(const std::string& property_name, unsigned int conformer_id = 0, SurfaceParams surface_params = SurfaceParams());
        
        /**
         * @brief Register a property factory for a given property name
         * 
         * This allows extensible addition of new properties without modifying core code.
         * 
         * @param property_name Name of the property (e.g., "esp", "hydrophobicity")
         * @param factory Function that creates the property object
         * 
         * @example
         *   molecule.registerProperty("esp", []() { return std::make_unique<SurfaceESPProperty>(); });
         */
        void registerProperty(const std::string& property_name, 
                            std::function<std::unique_ptr<ISurfaceProperty>()> factory);
        
        /**
         * @brief Get a computed property value by name and key
         * 
         * Generic method to retrieve any computed property without specialized getters.
         * 
         * @tparam T Type of the property value to retrieve
         * @param property_name Name of the property (e.g., "esp")
         * @param key Specific key within the property (e.g., "esp_charges", "electrostatic_potential")
         * @param conformer_id ID of the conformer (default: 0)
         * @return Reference to the property value
         * 
         * @throws std::runtime_error if property not computed or not found
         * 
         * @example
         *   auto charges = molecule.getProperty<std::vector<double>>("esp", "esp_charges");
         *   auto potentials = molecule.getProperty<std::vector<double>>("esp", "electrostatic_potential");
         */
        template<typename T>
        const T& getProperty(const std::string& property_name, const std::string& key, unsigned int conformer_id = 0) const {
            // Find the surface_id for this conformer (assuming it's the same as conformer_id for now)
            std::string cache_key = std::to_string(conformer_id) + "_" + property_name;
            auto cache_it = surface_property_cache_.find(cache_key);
            
            if (cache_it == surface_property_cache_.end()) {
                throw std::runtime_error("Property '" + property_name + "' not computed for conformer " + 
                                       std::to_string(conformer_id) + ". Call compute(\"" + property_name + "\") first.");
            }
            
            // Try to find the specific key in the cache
            auto key_it = cache_it->second.find(key);
            if (key_it == cache_it->second.end()) {
                throw std::runtime_error("Property key '" + key + "' not found for property '" + property_name + 
                                       "' in conformer " + std::to_string(conformer_id));
            }
            
            try {
                return std::any_cast<const T&>(key_it->second);
            } catch (const std::bad_any_cast& e) {
                throw std::runtime_error("Property key '" + key + "' has wrong type for property '" + property_name + "'");
            }
        }
        
        /**
         * @brief Register all available properties
         * 
         * Convenience method that registers all currently available properties.
         * This is useful for ensuring all properties are available without
         * manually registering each one.
         */
        void registerAllProperties();
        
        /**
         * @brief Get list of registered property names
         * 
         * @return Vector of strings containing all registered property names
         */
        std::vector<std::string> getRegisteredProperties() const;
        
        
        
        /**
         * @brief Get the number of vertices in a surface
         * 
         * @param conformer_id ID of the conformer (surface key)
         * @return Number of vertices
         * 
         * @throws std::runtime_error if surface is not found
         */
        unsigned int getSurfaceVertexCount(unsigned int conformer_id) const;
        
        /**
         * @brief Create or get a surface with the specified parameters
         *
         * This method creates a surface if it doesn't exist, or returns the existing one.
         * The surface is mapped to a specific conformer ID.
         *
         * @param conformer_id ID of the conformer to generate surface for (default: 0)
         * @param num_vertices Number of vertices to generate on the surface
         * @param radius Probe radius for surface generation in Angstroms (default: 1.2)
         * @param density Surface density parameter (default: 3.0)
         * @param hdensity High density parameter (default: 3.0)
         * @param type Surface type (default: "tses")
         * @param program Surface generation program (default: "msms")
         * @return Conformer ID used as the surface key
         *
         * @throws std::invalid_argument if parameters are invalid
         * @throws std::runtime_error if surface generation fails or conformer doesn't exist
         */
        int createSurface(unsigned int conformer_id = 0, unsigned int num_vertices = 1000, double radius = 1.2,
                         double density = 3.0, double hdensity = 3.0,
                         std::string type = "tses", std::string program = "msms", std::string sample_method = "full", bool addH = true);
        
        /**
         * @brief Get surface vertices for alignment
         * 
         * @param conformer_id ID of the conformer (surface key)
         * @return Reference to the surface vertices map
         * 
         * @throws std::runtime_error if surface is not found
         */
        const std::map<int, std::pair<std::string, std::vector<double>>>& getSurfaceVertices(unsigned int conformer_id) const;
        
        /**
         * @brief Clear cached surfaces
         * 
         * Removes all cached surface objects, forcing regeneration on next access.
         */
        void clearSurfaceCache();


        
        /**
         * @brief Apply a transformation to a conformer's coordinates
         * 
         * @param conformer_id ID of the conformer to transform
         * @param transformation 4x4 transformation matrix to apply
         * @return True if transformation was applied successfully
         * 
         * @throws std::runtime_error if conformer doesn't exist
         */
        bool transformConformer(unsigned int conformer_id, const Eigen::Isometry3d& transformation);
        
        /**
         * @brief Apply a transformation directly to surface vertices (efficient)
         * 
         * This method transforms the surface vertices directly without regenerating
         * the surface from scratch. Much more efficient than transformConformer + regenerate.
         * 
         * @param conformer_id ID of the conformer whose surface to transform
         * @param transformation 4x4 transformation matrix to apply
         * @return True if transformation was applied successfully
         * 
         * @throws std::runtime_error if conformer doesn't exist or surface not found
         */
        bool transformSurface(unsigned int conformer_id, const Eigen::Isometry3d& transformation);
        
        /**
         * @brief Align this molecule to another molecule using colored point cloud registration
         * 
         * Uses Open3D's colored ICP algorithm with surface properties as colors for better alignment.
         * The method computes surface properties (ESP) and uses them to color the point clouds
         * for more accurate registration.
         * 
         * @param target_mol The target molecule to align to
         * @param target_conformer_id Conformer ID of the target molecule's surface
         * @param source_conformer_id Conformer ID of this molecule's surface
         * @return 4x4 transformation matrix (Isometry3d) that aligns source to target
         * 
         * @throws std::runtime_error if conformers don't exist or surfaces not found
         */
        Eigen::Isometry3d align_to_molecule(Molecule& target_mol, unsigned int target_conformer_id = 0, unsigned int source_conformer_id = 0);
        
        /**
         * @brief Alignment result containing transformation and scores
         */
        struct AlignmentResult {
            Eigen::Isometry3d transformation;
            double fitness;        // ICP fitness score (lower is better)
            double inlier_rmse;    // Inlier RMSE (lower is better)
            bool used_colored_icp; // Whether colored ICP was used
        };
        
        /**
         * @brief Align this molecule to another molecule and return detailed alignment results
         * 
         * @param target_mol The target molecule to align to
         * @param target_conformer_id Conformer ID of the target molecule's surface
         * @param source_conformer_id Conformer ID of this molecule's surface
         * @return AlignmentResult containing transformation and scores
         */
        AlignmentResult align_to_molecule_detailed(Molecule& target_mol, unsigned int target_conformer_id = 0, unsigned int source_conformer_id = 0);
        
        /**
         * @brief Get the number of conformers available
         * 
         * @return Number of conformers
         */
        int getNumConformers() const;
        
        /**
         * @brief Get conformer IDs
         * 
         * @return Vector of available conformer IDs
         */
        std::vector<int> getConformerIds() const;
        
        /**
         * @brief Save molecule conformer to SDF format
         * 
         * @param conformer_id ID of the conformer to save
         * @param filename Output filename (should end with .sdf)
         * @param include_properties Whether to include molecular properties in the output
         * 
         * @throws std::runtime_error if conformer doesn't exist or file cannot be written
         */
        void saveConformerSDF(unsigned int conformer_id, const std::string& filename, bool include_properties = true) const;
        
        /**
         * @brief Save multiple conformers to SDF format with alignment scores
         * 
         * @param conformer_results Vector of conformer IDs and their alignment results
         * @param filename Output filename (should end with .sdf)
         * @param include_properties Whether to include molecular properties in the output
         * 
         * @throws std::runtime_error if conformers don't exist or file cannot be written
         */
        void saveMultiConformerSDF(const std::vector<std::pair<int, AlignmentResult>>& conformer_results, 
                                  const std::string& filename, bool include_properties = true) const;
        
        /**
         * @brief Save surface point cloud through molecule interface
         * 
         * @param conformer_id ID of the conformer whose surface to save
         * @param filename Output filename (should end with .ply)
         * @param include_normals Whether to include normal vectors
         * @param include_colors Whether to include color information
         * 
         * @throws std::runtime_error if surface doesn't exist or file cannot be written
         */
        void saveSurfacePointCloud(unsigned int conformer_id, const std::string& filename, bool include_normals = true, bool include_colors = false, const std::vector<std::string>& properties = {}) const;
        
        /**
         * @brief Save surface mesh through molecule interface
         * 
         * @param conformer_id ID of the conformer whose surface to save
         * @param filename Output filename (should end with .ply)
         * @param include_normals Whether to include normal vectors
         * @param include_colors Whether to include color information
         * 
         * @throws std::runtime_error if surface doesn't exist or file cannot be written
         */
        void saveSurfaceMesh(unsigned int conformer_id, const std::string& filename, bool include_normals = true, bool include_colors = false, const std::vector<std::string>& properties = {}) const;
        
        /**
         * @brief Load molecules from a multi-molecule file (SDF)
         * 
         * @param filename Path to the SDF file
         * @return Vector of unique_ptr to Molecule objects, one for each molecule in the file
         * 
         * @throws std::runtime_error if file cannot be read or parsed
         */
        static std::vector<std::unique_ptr<Molecule>> loadFromMultiMolFile(const std::string& filename);
        
        /**
         * @brief Get file extension from filename
         * 
         * @param filename Input filename
         * @return Lowercase file extension (e.g., "sdf", "mol", "pdb")
         */
        static std::string getFileExtension(const std::string& filename);
        
        /**
         * @brief Add hydrogens to molecule only if none are present
         * 
         * Checks if hydrogens are already present in the molecule before adding them.
         * This prevents duplicate hydrogens when loading molecules from files that
         * already contain hydrogen atoms.
         */
        void add_hydrogens_to_mol_if_needed();
        
        /**
         * @brief Remove explicit hydrogens from the molecule
         * 
         * Uses RDKit's removeHs function to remove explicit hydrogen atoms from the
         * molecular structure.
         */
        void remove_hydrogens_from_mol();

        /**
         * @brief Set the name of the molecule
         */
        void set_name(const std::string& n) { name = n; }

        /**
         * @brief Get the name of the molecule
         */
        std::string get_name() const { return name; }

    private:
        /**
         * @brief Name of the molecule
         */
        std::string name;

        /**
         * @brief Original SMILES string
         */
        std::string smiles;
        
        /**
         * @brief Formal charge of the molecule
         */
        int charge = 0;
        
        /**
         * @brief RDKit molecule object (read-write)
         */
        std::unique_ptr<RDKit::RWMol> mol;
        
        /**
         * @brief Cache for generated surfaces
         *
         * Maps conformer IDs to cached Surface objects for efficient reuse.
         * Each conformer has its own surface representation.
         */
        mutable std::map<int, std::shared_ptr<Surface>> surface_cache_;
        
        /**
         * @brief Cache for computed surface properties
         *
         * Maps surface_key + property_name to property cache for efficient reuse.
         */
        mutable std::map<std::string, std::unordered_map<std::string, std::any>> surface_property_cache_;
        
        /**
         * @brief Registry of property factories
         *
         * Maps property names to factory functions that create property objects.
         */
        std::map<std::string, std::function<std::unique_ptr<ISurfaceProperty>()>> property_registry_;
        
        /**
         * @brief Register default properties (ESP, etc.)
         * 
         * This method registers commonly used properties so they're available
         * without manual registration.
         */
        void registerDefaultProperties();
        

        /**
         * @brief Parse SMILES string and create RDKit molecule object
         * 
         * Converts the SMILES string to an RDKit RWMol object using RDKit's
         * SMILES parser. Sets mol to nullptr if parsing fails.
         */
        void save_smiles_to_mol();
        
        /**
         * @brief Calculate and store formal charge from RDKit molecule
         * 
         * Uses RDKit's MolOps::getFormalCharge to compute the formal charge
         * of the molecule and stores it in the charge member variable.
         */
        void get_charge_from_mol();
        
        /**
         * @brief Generate 3D coordinates for the molecule
         * 
         * Uses RDKit's standard conformer generation with force field optimization
         * if the molecule doesn't already have valid 3D structure.
         */
        void generate3DCoordinates();
        
        /**
         * @brief Generate a single optimized conformer
         * 
         * Uses RDKit + force field optimization to generate one high-quality conformer.
         * This is useful when you only need one conformer instead of multiple.
         */
        void generateSingleConformer();
        
        /**
         * @brief Add explicit hydrogens to the molecule
         * 
         * Uses RDKit's addHs function to add explicit hydrogen atoms to the
         * molecular structure. Currently disabled for surface generation.
         */
        void add_hydrogens_to_mol();
        
        
        /**
         * @brief Generate a cache key for surface parameters
         * 
         * @param num_vertices Number of vertices
         * @param radius Probe radius
         * @param density Surface density
         * @param hdensity High density
         * @param type Surface type
         * @param program Surface program
         * @return String key for surface cache
         */
        std::string generateSurfaceCacheKey(unsigned int num_vertices, double radius, double density, 
                                          double hdensity, const std::string& type, 
                                          const std::string& program) const;
};