/**
 * @file surface.h
 * @brief Surface generation and manipulation for molecular structures
 * @author Sebastian
 * @date 2025
 */

#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <unordered_map>
#include <any>
#include <memory>
#include <Eigen/Geometry>

#include <GraphMol/GraphMol.h>

#include "molecule.h"
#include "props/surface_property.h"
#include "surface_fibonacci.h"

#ifndef SURFACE_H
#define SURFACE_H

/**
 * @class Surface
 * @brief Molecular surface representation and property computation
 * 
 * The Surface class represents a triangulated molecular surface generated from a
 * molecule's 3D coordinates. It provides functionality for surface property computation,
 * vertex sampling, and integration with various surface generation programs.
 * 
 * Key features:
 * - Surface representation with vertices and faces
 * - Multiple surface generation program support (MSMS, etc.)
 * - Van der Waals radius database with 118 elements
 * - Surface property computation and caching
 * - Vertex sampling algorithms (random, FPS)
 * 
 * @see Molecule
 * @see ISurfaceProperty
 */
class Surface {
public:
    /**
     * @brief Construct a new Surface object
     * 
     * @param molecule Reference to the Molecule object containing the molecular structure
     * @param conformer_id ID of the conformer to generate surface for
     * @param num_vertices Number of vertices to generate on the surface (must be > 4)
     * @param radius Probe radius for surface generation in Angstroms (default: 1.2)
     * @param density Surface density parameter (default: 3.0)
     * @param hdensity High density parameter (default: 3.0)
     * @param type Surface type - "tses" for triangulated solvent excluded surface, "ases" for analytical (default: "tses")
     * @param program Surface generation program to use - "msms" or "fibonacci" (default: "msms")
     * @param sample_method Vertex sampling method - "random", "fps" (farthest point sampling), or "full" (default: "full")
     * 
     * @throws std::invalid_argument if num_vertices <= 4, radius <= 0, density <= 0, or hdensity <= 0
     * @throws std::invalid_argument if type is not "tses" or "ases"
     * @throws std::invalid_argument if sample_method is not "random", "fps", or "full"
     * @throws std::runtime_error if conformer doesn't exist
     */
    Surface(const Molecule& molecule, int conformer_id, int num_vertices, double radius = 1.2, double density = 3.0, double hdensity = 3.0, std::string type = "tses", std::string program = "msms", std::string sample_method = "full");
    
    /**
     * @brief Destroy the Surface object
     */
    ~Surface();

    /**
     * @brief Main surface generation method - handles all input types automatically
     * 
     * This method automatically detects the input file type and handles the complete pipeline:
     * - If input has .xyzr extension: directly uses the xyzr file
     * - Otherwise: generates 3D coordinates from molecule, creates xyzr file, then generates surface
     * 
     * @param input_file Input file path (with or without extension)
     * @param output_file Output file prefix for generated surface files
     * 
     * @throws std::runtime_error if surface generation fails
     * @throws std::runtime_error if molecule is invalid or has no 3D coordinates
     */
    void generateSurface(std::string_view input_file, std::string_view output_file);
    
    /**
     * @brief Legacy method for direct xyzr file surface generation
     * 
     * Directly generates surface from an existing xyzr file without any preprocessing.
     * The xyzr file must be in the standard format: element x y z radius
     * 
     * @param xyzr_file Path to the xyzr file
     * @param output_file Output file prefix for generated surface files
     * 
     * @throws std::runtime_error if xyzr file cannot be opened or is malformed
     * @throws std::runtime_error if surface generation program fails
     */
    void generateSurfaceFromXYZR(std::string_view xyzr_file, std::string_view output_file);
    
    /**
     * @brief Generate surface using Fibonacci sphere sampling method
     * 
     * This method uses the Fibonacci sphere sampling approach with Open3D's ball pivoting
     * algorithm to generate a molecular surface. It's independent of external programs
     * like MSMS and provides an alternative surface generation approach.
     * 
     * @param fibonacci_params Parameters for Fibonacci surface generation
     * 
     * @throws std::runtime_error if molecule has no conformers or surface generation fails
     */
    void generateFibonacciSurface(const spalt::FibonacciSurfaceParams& fibonacci_params = spalt::FibonacciSurfaceParams());
    
    /**
     * @brief Select a subset of vertices using the specified sampling method
     * 
     * @param n_vertices Number of vertices to select
     * @param sampling_method Sampling algorithm: "random", "fps" (farthest point sampling)
     * 
     * @return Map of selected vertices with their indices, atom names, and coordinates/normals
     * 
     * @throws std::invalid_argument if sampling_method is not recognized
     */
    std::map<int, std::pair<std::string, std::vector<double>>> select_n_vertices(size_t n_vertices, std::string sampling_method);

    /**
     * @brief Get the number of vertices in the surface
     * 
     * @return Number of vertices
     */
    int get_num_vertices() const;

    /**
     * @brief Register a surface property for computation
     * 
     * @param prop Unique pointer to the surface property object
     */
    void registerProperty(std::unique_ptr<ISurfaceProperty> prop);
    
    /**
     * @brief Compute a specific surface property
     * 
     * @param key Key identifying the property to compute
     * 
     * @throws std::runtime_error if property is not registered
     */
    void computeProperty(const std::string& key);
    
    /**
     * @brief Get a computed surface property by key
     * 
     * @tparam T Type of the property value
     * @param key Key identifying the property
     * @return Reference to the property value
     * 
     * @throws std::runtime_error if property has not been computed
     * @throws std::bad_any_cast if property type doesn't match T
     */
    template<typename T> const T& getProperty(const std::string& key) const {
        auto it = cache_.find(key);
        if (it == cache_.end())
            throw std::runtime_error("Property not computed: " + key);
        return std::any_cast<const T&>(it->second);
    }
    
    /**
     * @brief Invalidate all cached surface properties
     * 
     * Clears the property cache, forcing recomputation of all properties on next access.
     */
    void invalidateProperties();
    
    /**
     * @brief Reference to the associated Molecule object
     */
    const Molecule& molecule;
    
    /**
     * @brief ID of the conformer this surface represents
     */
    int conformer_id;

    /**
     * @brief Surface vertices map: vertex_index -> (atom_name, {coordinates, normals})
     * 
     * Each vertex contains:
     * - Atom name or identifier
     * - Vector with 6 elements: [x, y, z, nx, ny, nz] (coordinates and normal)
     * 
     * Normal Handling Strategy:
     * - MSMS surfaces: Normals are computed by MSMS and stored in elements 3-5
     * - Fibonacci surfaces: Normals are computed from mesh faces using face averaging
     * - Alignment: Existing normals are preserved; if missing, Open3D estimates them
     * - Export: Normals are included in PLY files when requested
     */
    std::map<int, std::pair<std::string, std::vector<double>>> vertices = {};

    
    /**
     * @brief Apply a transformation directly to the surface vertices
     * 
     * This is much more efficient than regenerating the surface from scratch.
     * Transforms both vertex coordinates and normals.
     * 
     * @param transformation 4x4 transformation matrix to apply
     */
    void transformSurface(const Eigen::Isometry3d& transformation);
    
    /**
     * @brief Save surface vertices as a point cloud in PLY format
     * 
     * @param filename Output filename (should end with .ply)
     * @param include_normals Whether to include normal vectors in the output
     * @param include_colors Whether to include color information (default: false)
     * @param properties Vector of property names to include as color columns
     * @param property_cache Property cache containing computed values
     * 
     * @throws std::runtime_error if file cannot be written
     */
    void savePointCloud(const std::string& filename, bool include_normals = true, bool include_colors = false, const std::vector<std::string>& properties = {}, const std::unordered_map<std::string, std::any>* property_cache = nullptr) const;
    
    /**
     * @brief Save surface mesh in PLY format (vertices + faces)
     * 
     * @param filename Output filename (should end with .ply)
     * @param include_normals Whether to include normal vectors in the output
     * @param include_colors Whether to include color information (default: false)
     * @param properties Vector of property names to include as color columns
     * @param property_cache Property cache containing computed values
     * 
     * @throws std::runtime_error if file cannot be written
     */
    void saveMesh(const std::string& filename, bool include_normals = true, bool include_colors = false, const std::vector<std::string>& properties = {}, const std::unordered_map<std::string, std::any>* property_cache = nullptr) const;
private:
    /**
     * @brief Convert property value to RGB color
     * 
     * @param property_name Name of the property
     * @param vertex_index Index of the vertex
     * @param fallback_atom_name Atom name to use for fallback coloring
     * @return RGB color as tuple (red, green, blue) in range [0, 255]
     */
    std::tuple<int, int, int> getVertexColor(const std::string& property_name, int vertex_index, const std::string& fallback_atom_name, const std::unordered_map<std::string, std::any>* property_cache) const;
    
    /**
     * @brief Get normalized color value for vertex (-1 to 1)
     * 
     * @param property_name Name of the property
     * @param vertex_index Index of the vertex
     * @param fallback_atom_name Atom name to use for fallback coloring
     * @param property_cache Property cache containing computed values
     * @return Normalized color value in range [-1, 1]
     */
    double getVertexColorValue(const std::string& property_name, int vertex_index, const std::string& fallback_atom_name, const std::unordered_map<std::string, std::any>* property_cache) const;
    
    /**
     * @brief Normalize property value to [-1, 1] range
     * 
     * @param property_name Name of the property
     * @param value Property value to normalize
     * @param property_cache Property cache containing min/max values
     * @return Normalized value in range [-1, 1]
     */
    double normalizePropertyValue(const std::string& property_name, double value, const std::unordered_map<std::string, std::any>* property_cache) const;
    
    /**
     * @brief Get color for atom type (CPK coloring scheme)
     * 
     * @param atom_name Name of the atom
     * @return RGB color as tuple (red, green, blue) in range [0, 255]
     */
    std::tuple<int, int, int> getAtomTypeColor(const std::string& atom_name) const;
    
    /**
     * @brief Convert property value to RGB color
     * 
     * @param property_name Name of the property
     * @param value Property value
     * @return RGB color as tuple (red, green, blue) in range [0, 255]
     */
    std::tuple<int, int, int> propertyValueToColor(const std::string& property_name, double value, const std::unordered_map<std::string, std::any>* property_cache) const;
    
    /**
     * @brief Van der Waals radii database (Grimme) in Angstroms
     * 
     * Contains van der Waals radii for 118 elements based on Grimme's parameterization.
     * Used for automatic xyzr file generation from 3D coordinates.
     */
    static const std::map<std::string, double> vdw_radii_;

    /**
     * @brief Number of vertices to generate on the surface
     */
    int num_vertices;
    
    /**
     * @brief Probe radius for surface generation in Angstroms
     */
    double radius;
    
    /**
     * @brief Surface density parameter for triangulation
     */
    double density;
    
    /**
     * @brief High density parameter for surface generation
     */
    double hdensity;
    
    /**
     * @brief Surface type - "tses" for triangulated, "ases" for analytical
     */
    std::string type;
    
    /**
     * @brief Surface generation program to use (e.g., "msms")
     */
    std::string program;
    
    /**
     * @brief Minimum number of columns expected in xyzr files
     */
    static constexpr size_t EXPECTED_MIN_COLUMNS = 4;
    
    /**
     * @brief Expected number of columns in vertex files from surface generation programs
     */
    static constexpr size_t EXPECTED_VERT_COLUMNS = 10;
    
    /**
     * @brief Default van der Waals radius for unknown elements
     */
    static constexpr size_t DEFAULT_VDW_RADIUS = 1.7;
    
    /**
     * @brief Surface faces as triangles (vertex indices)
     */
    std::vector<std::vector<int>> faces = {};

    /**
     * @brief Registered surface properties
     */
    std::unordered_map<std::string, std::unique_ptr<ISurfaceProperty>> props_;
    
    /**
     * @brief Cache for computed surface property values
     */
    std::unordered_map<std::string, std::any> cache_;
    
    /**
     * @brief Compute vertex normal from mesh faces
     * 
     * @param vertex_index Index of the vertex
     * @param vertices Vector of vertex positions
     * @param faces Vector of face indices
     * @return Computed vertex normal as unit vector
     */
    Eigen::Vector3d computeVertexNormal(size_t vertex_index, 
                                      const std::vector<Eigen::Vector3d>& vertices, 
                                      const std::vector<Eigen::Vector3i>& faces) const;

    /**
     * @brief Run MSMS surface generation program
     * 
     * @param xyzr_file Input xyzr file path
     * @param output_file Output file prefix
     * 
     * @throws std::runtime_error if MSMS command fails
     */
    void run_msms(std::string_view xyzr_file, std::string_view output_file);

    /**
     * @brief Read vertex file generated by surface generation program
     * 
     * @param vert_file Path to vertex file
     * @return true if successful, false otherwise
     */
    bool read_vert_file(std::string_view vert_file);

    /**
     * @brief Read face file generated by surface generation program
     * 
     * @param face_file Path to face file
     * @return true if successful, false otherwise
     */
    bool read_face_file(std::string_view face_file);

    /**
     * @brief Select vertices using random sampling
     * 
     * @param n_vertices Number of vertices to select
     * @return Map of selected vertices
     */
    std::map<int, std::pair<std::string, std::vector<double>>> select_n_vertices_random(size_t n_vertices);

    /**
     * @brief Select vertices using farthest point sampling (FPS)
     * 
     * @param n_vertices Number of vertices to select
     * @return Map of selected vertices
     */
    std::map<int, std::pair<std::string, std::vector<double>>> select_n_vertices_fps(size_t n_vertices);
    
    /**
     * @brief Generate 3D coordinates for the molecule if not already present
     * 
     * Uses RDKit's conformer generation to create 3D coordinates from 2D structure.
     * Only generates coordinates if the molecule doesn't already have valid 3D coordinates.
     * 
     * @throws std::runtime_error if molecule is invalid or conformer generation fails
     */
    void generate3DCoordinates();
    
    /**
     * @brief Generate xyzr file from molecule's 3D coordinates
     * 
     * Creates an xyzr file in the format: element x y z radius
     * Uses the van der Waals radii database for proper radius assignment.
     * 
     * @param output_file Output xyzr file path
     * 
     * @throws std::runtime_error if molecule has no conformers or file creation fails
     */
    void generateXYZRFile(std::string_view output_file);
    
    /**
     * @brief Get file extension from filename
     * 
     * @param filename Input filename
     * @return File extension (without dot) or empty string if no extension
     */
    std::string getFileExtension(std::string_view filename) const;
    
    /**
     * @brief Apply vertex subsampling based on the sample_method
     * 
     * This method applies the specified sampling method to reduce the number
     * of vertices to the target num_vertices. If sample_method is "full",
     * no subsampling is applied.
     */
    void applySubsampling();

private:
    // Member variables
    std::string sample_method;          ///< Vertex sampling method (random/fps/full)
    
public:
};

#endif