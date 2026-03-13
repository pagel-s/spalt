#ifndef CLI_H
#define CLI_H

#include <filesystem>
#include <string>
#include <vector>
#include "conformer_params.h"
#include "props/esp.h"

/**
 * @brief Command Line Interface for molecular surface alignment
 *
 * Supports aligning molecules from various sources (files, SMILES, folders)
 * to a reference molecule and exporting aligned surfaces/conformers.
 */
class CLI {
  public:
    /**
     * @brief Parse command line arguments and execute alignment
     *
     * @param argc Number of command line arguments
     * @param argv Array of command line arguments
     * @return Exit code (0 for success, non-zero for error)
     */
    static int run(int argc, char* argv[]);

    /**
     * @brief Print usage information
     *
     * @param program_name Name of the executable
     */
    static void printUsage(const std::string& program_name);

    /**
     * @brief Print version information
     */
    static void printVersion();

    /**
     * @brief Parse surface parameters from command line arguments
     *
     * @param argc Number of command line arguments
     * @param argv Array of command line argument strings
     * @param start_index Index to start parsing from (will be updated to point after parsed
     * arguments)
     * @return SurfaceParams struct with parsed parameters and defaults for unspecified options
     */
    static SurfaceParams parseSurfaceParams(int argc, char* argv[], int& start_index);

    /**
     * @brief Surface generation parameters for CLI operations
     *
     * Parameter structure for command-line interface operations.
     * Contains surface generation and conformer generation parameters with
     * sensible defaults for typical molecular alignment workflows.
     */
    struct SurfaceParams {
        int num_vertices =
            1000;  ///< Number of vertices to generate on the molecular surface (default: 1000)
        double radius = 1.2;    ///< Probe radius for surface generation in Angstroms (default: 1.2)
        double density = 3.0;   ///< Surface density parameter for MSMS algorithm (default: 3.0)
        double hdensity = 3.0;  ///< High density parameter for MSMS algorithm (default: 3.0)
        std::string type = "tses";  ///< Surface type: "tses" (triangulated solvent excluded
                                    ///< surface) or "ases" (analytical surface) (default: "tses")
        std::string program = "fibonacci";  ///< Surface generation program: "msms" or "fibonacci"
                                            ///< (default: "fibonacci")
        std::string sample_method =
            "full";  ///< Vertex sampling method: "full" (all vertices), "fps" (farthest point
                     ///< sampling), or "random" (default: "full")
        std::vector<std::string> properties = {
            "esp", "hydrophobicity", "hbond"};  ///< List of surface properties to compute (default:
                                                ///< compute all available properties)
        std::string output_type = "mesh";  ///< Output format: "mesh" (triangular mesh) or "points"
                                           ///< (point cloud) (default: "mesh")
        bool include_normals = false;  ///< Whether to include vertex normals in output (default:
                                       ///< false for better PyMOL compatibility)
        int top_n = -1;  ///< Number of top conformers to save (-1 = save all, >0 = save only top N
                         ///< conformers) (default: -1)

        // Advanced conformer generation parameters
        ConformerParams conformer_params;  ///< Advanced conformer generation settings including
                                           ///< clustering and optimization parameters
        bool save_meshes =
            false;  ///< Whether to save individual mesh files for each conformer (default: false)

        // ESP charge calculation method
        ChargeMethod charge_method =
            ChargeMethod::RDKIT;  ///< Charge calculation method for ESP: XTB (quantum chemistry) or
                                  ///< RDKIT (Gasteiger charges) (default: RDKIT)

        // Hydrogen control parameters
        bool addH = true;  ///< Whether to add explicit hydrogens to molecules (default: true)
    };

  private:
    /**
     * @brief Align molecules from a file to a reference molecule
     *
     * @param reference_file Path to reference molecule file
     * @param input_file Path to input molecule file (can be multi-molecule SDF)
     * @param output_dir Output directory for aligned files
     * @param surface_params Surface generation parameters
     */
    static void alignFromFile(const std::string& reference_file, const std::string& input_file,
                              const std::string& output_dir, const SurfaceParams& surface_params);

    /**
     * @brief Align molecules from SMILES to a reference molecule
     *
     * @param reference_file Path to reference molecule file
     * @param smiles SMILES string to align
     * @param output_dir Output directory for aligned files
     * @param surface_params Surface generation parameters
     */
    static void alignFromSMILES(const std::string& reference_file, const std::string& smiles,
                                const std::string& output_dir, const SurfaceParams& surface_params);

    /**
     * @brief Align molecules from a folder to a reference molecule
     *
     * @param reference_file Path to reference molecule file
     * @param input_dir Input directory containing molecule files
     * @param output_dir Output directory for aligned files
     * @param surface_params Surface generation parameters
     */
    static void alignFromFolder(const std::string& reference_file, const std::string& input_dir,
                                const std::string& output_dir, const SurfaceParams& surface_params);

    /**
     * @brief Create output directory if it doesn't exist
     *
     * @param output_dir Path to output directory (will be created recursively if needed)
     */
    static void createOutputDir(const std::string& output_dir);

    /**
     * @brief Get all molecule files from a directory
     *
     * @param dir_path Directory path to search for molecule files
     * @return Vector of molecule file paths (supports .sdf, .mol, .mol2, .pdb extensions)
     */
    static std::vector<std::string> getMoleculeFiles(const std::string& dir_path);
};

#endif  // CLI_H
