#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "molecule.h"
#include "props/esp.h"
#include "surface.h"

void printUsage(const char* program_name) {
    std::cout << "SPALT - Molecular Surface Alignment Tool\n" << std::endl;
    std::cout << "Usage: " << program_name << " <reference> <input> <output_dir> [options]\n"
              << std::endl;

    std::cout << "Arguments:" << std::endl;
    std::cout << "  reference         Reference molecule - one of:" << std::endl;
    std::cout << "                      - Path to SDF file" << std::endl;
    std::cout << "                      - SMILES string (will generate 3D coordinates)"
              << std::endl;
    std::cout << "  input             Input molecule(s) - one of:" << std::endl;
    std::cout << "                      - Path to SDF file (single or multi-molecule)" << std::endl;
    std::cout << "                      - Path to directory containing SDF files" << std::endl;
    std::cout << "                      - Path to text file containing SMILES (one per line)"
              << std::endl;
    std::cout << "                      - SMILES string (will generate 3D coordinates)"
              << std::endl;
    std::cout << "  output_dir        Directory to save aligned surfaces and conformers\n"
              << std::endl;

    std::cout << "Options:" << std::endl;
    std::cout << "  --vertices N      Number of surface vertices (default: 1000)" << std::endl;
    std::cout << "  --radius R        Probe radius in Angstroms (default: 1.2)" << std::endl;
    std::cout << "  --density D       Surface density (default: 3.0)" << std::endl;
    std::cout << "  --hdensity H      High density (default: 3.0)" << std::endl;
    std::cout << "  --type T          Surface type: tses or ases (default: tses)" << std::endl;
    std::cout << "  --mesh M          Surface method: msms, fibonacci (default: fibonacci)" << std::endl;
    std::cout << "  --sample-method S Vertex sampling method: full, fps, random (default: full)" << std::endl;
    std::cout << "  --top-n N          Save only the top N best aligned conformers (default: all)" << std::endl;
    std::cout
        << "  --properties P    Surface properties to compute (default: esp,hydrophobicity,hbond,none):"
        << std::endl;
    std::cout << "                      - esp: Electrostatic potential" << std::endl;
    std::cout << "                      - hb: Hydrogen bond potential (abbreviation for hbond)" << std::endl;
    std::cout << "                      - hy: Hydrophobicity values (abbreviation for hydrophobicity)" << std::endl;
    std::cout << "                      - hbond: Hydrogen bond potential (full name)" << std::endl;
    std::cout << "                      - hydrophobicity: Hydrophobicity values (full name)" << std::endl;
    std::cout << "                      - all: Compute all available properties" << std::endl;
    std::cout << "                      - Comma-separated list: esp,hydrophobicity" << std::endl;
    std::cout << "                      - none: Skip property computation" << std::endl;
    std::cout << "  --output-type T   Surface output format (default: mesh):" << std::endl;
    std::cout << "                      - mesh: Triangulated mesh (recommended for PyMOL)"
              << std::endl;
    std::cout << "                      - points: Point cloud" << std::endl;
    std::cout << "  --include-normals Include normal vectors in output (default: false)"
              << std::endl;
    std::cout
        << "  --conformers N    [DEPRECATED] Number of diverse conformers to select (maps to --num-clusters)"
        << std::endl;
    std::cout << "                      Use --num-clusters instead for clarity (default: 3)"
              << std::endl;
    std::cout
        << "  --save-meshes     Save individual mesh files for each conformer (default: false)"
        << std::endl;
    std::cout << "  --total-confs N   Total conformers to generate before clustering (default: 50)"
              << std::endl;
    std::cout << "  --num-clusters N  Number of diverse conformers to select via clustering (default: 5)"
              << std::endl;
    std::cout
        << "  --use-advanced    Use advanced conformer generation with clustering (default: false)"
        << std::endl;
    std::cout
        << "  --random-seed N   Random seed for reproducible conformer generation (default: 4)"
        << std::endl;
    std::cout
        << "  --charge-method M Charge calculation method for ESP (default: rdkit):"
        << std::endl;
    std::cout << "                      - rdkit: Fast Gasteiger charges (built-in)"
              << std::endl;
    std::cout << "                      - xtb: Quantum chemistry charges (requires xtb program)"
              << std::endl;
    std::cout << "  --addH            Add explicit hydrogens to molecules (default: true)"
              << std::endl;
    std::cout << "  --removeH         Remove explicit hydrogens from molecules (default: false)\n"
              << std::endl;

    std::cout << "Examples:" << std::endl;
    std::cout << "  " << program_name << " reference.sdf input.sdf ./output/" << std::endl;
    std::cout << "  " << program_name << " reference.sdf ./molecules/ ./output/" << std::endl;
    std::cout << "  " << program_name
              << " \"C1CCCCC1\" \"O=C1NC2=C(C1=CC3=C(CCNC4=O)C4=CN3)C=C(Br)C=C2\" ./output/"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf multi_molecules.sdf ./output/ --vertices 2000" << std::endl;
    std::cout << "  " << program_name << " reference.sdf input.sdf ./output/ --properties all"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf input.sdf ./output/ --properties esp,hy,hb" << std::endl;
    std::cout << "  " << program_name << " reference.sdf input.sdf ./output/ --properties none"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf input.sdf ./output/ --output-type mesh --include-normals"
              << std::endl;
    std::cout << "  " << program_name << " reference.sdf input.sdf ./output/ --output-type points"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf \"CCO\" ./output/ --conformers 5 --save-meshes" << std::endl;
    std::cout << "  " << program_name << " reference.sdf smiles.txt ./output/ --conformers 2"
              << std::endl;
    std::cout
        << "  " << program_name
        << " reference.sdf \"CCO\" ./output/ --use-advanced --total-confs 100 --num-clusters 8"
        << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf \"CCO\" ./output/ --use-advanced --random-seed 42 --save-meshes"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf \"CCO\" ./output/ --charge-method xtb --properties esp"
              << std::endl;
    std::cout << "  " << program_name
              << " reference.sdf \"CCO\" ./output/ --removeH --properties hb"
              << std::endl;
}

std::vector<std::string> parseProperties(const std::string& properties_str) {
    std::vector<std::string> properties;

    if (properties_str == "all") {
        properties = {"esp", "hydrophobicity", "hbond"};
    } else if (properties_str == "none") {
        properties = {};  // Empty vector means no properties
    } else {
        // Parse comma-separated list
        std::stringstream ss(properties_str);
        std::string property;
        while (std::getline(ss, property, ',')) {
            // Trim whitespace
            property.erase(0, property.find_first_not_of(" \t"));
            property.erase(property.find_last_not_of(" \t") + 1);

            if (!property.empty()) {
                // Handle abbreviated property names
                if (property == "esp") {
                    properties.push_back("esp");
                } else if (property == "hb") {
                    properties.push_back("hbond");
                } else if (property == "hy") {
                    properties.push_back("hydrophobicity");
                } else {
                    // Use the property name as-is (supports full names too)
                    properties.push_back(property);
                }
            }
        }
    }

    return properties;
}

SurfaceParams parseSurfaceParams(int argc, char* argv[], int& start_index) {
    SurfaceParams params;
    
    // Track parameter usage for conflict detection
    bool conformers_specified = false;
    bool num_clusters_specified = false;

    while (start_index < argc) {
        std::string arg = argv[start_index];

        if (arg == "--vertices" && start_index + 1 < argc) {
            params.num_vertices = std::stoi(argv[++start_index]);
        } else if (arg == "--radius" && start_index + 1 < argc) {
            params.radius = std::stod(argv[++start_index]);
        } else if (arg == "--density" && start_index + 1 < argc) {
            params.density = std::stod(argv[++start_index]);
        } else if (arg == "--hdensity" && start_index + 1 < argc) {
            params.hdensity = std::stod(argv[++start_index]);
        } else if (arg == "--type" && start_index + 1 < argc) {
            params.type = argv[++start_index];
        } else if (arg == "--mesh" && start_index + 1 < argc) {
            std::string mesh = argv[++start_index];
            if (mesh == "msms" || mesh == "m" || mesh == "ms") {
                params.program = "msms";
            } else if (mesh == "fib" || mesh == "f" || mesh == "fibonacci") {
                params.program = "fibonacci";
            } else {
                std::cerr << "Error: Invalid mesh method '" << mesh << "'. Must be 'msms'[, 'm', 'ms'], 'fibonacci'[, 'f', or 'fib']." << std::endl;
                return SurfaceParams{};
            }
        } else if (arg == "--sample-method" && start_index + 1 < argc) {
            std::string method = argv[++start_index];
            if (method == "full" || method == "fps" || method == "random") {
                params.sample_method = method;
            } else {
                std::cerr << "Error: Invalid sample method '" << method << "'. Must be 'full', 'fps', or 'random'." << std::endl;
                return SurfaceParams{};
            }
        } else if (arg == "--top-n" && start_index + 1 < argc) {
            int top_n = std::stoi(argv[++start_index]);
            if (top_n <= 0) {
                std::cerr << "Error: --top-n must be a positive integer." << std::endl;
                return SurfaceParams{};
            }
            params.top_n = top_n;
        } else if (arg == "--properties" && start_index + 1 < argc) {
            params.properties = parseProperties(argv[++start_index]);
        } else if (arg == "--output-type" && start_index + 1 < argc) {
            std::string output_type = argv[++start_index];
            if (output_type == "mesh" || output_type == "points") {
                params.output_type = output_type;
            } else {
                std::cerr << "Warning: Invalid output type '" << output_type << "'. Using 'mesh'."
                          << std::endl;
            }
        } else if (arg == "--include-normals") {
            params.include_normals = true;
        } else if (arg == "--conformers" && start_index + 1 < argc) {
            if (num_clusters_specified) {
                std::cerr << "Warning: --conformers specified after --num-clusters. "
                          << "Using --num-clusters value (" << params.conformer_params.num_clusters 
                          << ") and ignoring --conformers (" << argv[start_index + 1] << ")" << std::endl;
            } else {
                conformers_specified = true;
                params.conformer_params.num_clusters = std::stoi(argv[++start_index]);
                std::cerr << "Info: --conformers maps to num_clusters for backward compatibility. "
                          << "Consider using --num-clusters instead for clarity." << std::endl;
            }
            ++start_index; // Skip the value even if we ignored it
        } else if (arg == "--save-meshes") {
            params.save_meshes = true;
        } else if (arg == "--total-confs" && start_index + 1 < argc) {
            params.conformer_params.total_conformers = std::stoi(argv[++start_index]);
        } else if (arg == "--num-clusters" && start_index + 1 < argc) {
            if (conformers_specified) {
                std::cerr << "Warning: --num-clusters specified after --conformers. "
                          << "Overriding --conformers value with --num-clusters (" << argv[start_index + 1] << ")" << std::endl;
            }
            num_clusters_specified = true;
            params.conformer_params.num_clusters = std::stoi(argv[++start_index]);
        } else if (arg == "--use-advanced") {
            params.conformer_params.enable_optimization = true;
        } else if (arg == "--random-seed" && start_index + 1 < argc) {
            params.conformer_params.random_seed = std::stoi(argv[++start_index]);
        } else if (arg == "--charge-method" && start_index + 1 < argc) {
            std::string method = argv[++start_index];
            if (method == "xtb") {
                params.charge_method = ChargeMethod::XTB;
            } else if (method == "rdkit") {
                params.charge_method = ChargeMethod::RDKIT;
            } else {
                std::cerr << "Warning: Invalid charge method '" << method 
                          << "'. Must be 'xtb' or 'rdkit'. Using default 'rdkit'." << std::endl;
                params.charge_method = ChargeMethod::RDKIT;
            }
        } else if (arg == "--addH") {
            params.addH = true;
        } else if (arg == "--removeH") {
            params.addH = false;
        } else {
            break;  // Unknown argument, stop parsing
        }

        ++start_index;
    }
    
    // Final validation and info
    if (params.conformer_params.enable_optimization) {
        std::cerr << "Info: Advanced conformer generation enabled:" << std::endl;
        std::cerr << "  - Total conformers to generate: " << params.conformer_params.total_conformers << std::endl;
        std::cerr << "  - Clusters for diverse selection: " << params.conformer_params.num_clusters << std::endl;
        std::cerr << "  - Final result: up to " << params.conformer_params.num_clusters 
                  << " diverse conformers selected from " << params.conformer_params.total_conformers << " generated" << std::endl;
    }

    return params;
}

std::vector<std::string> loadSmilesFromFile(const std::string& file_path) {
    std::vector<std::string> smiles_list;
    std::ifstream file(file_path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open SMILES file: " + file_path);
    }

    std::string line;
    int line_number = 0;
    while (std::getline(file, line)) {
        line_number++;

        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r"));
        line.erase(line.find_last_not_of(" \t\r") + 1);

        // Skip empty lines and comments
        if (!line.empty() && line[0] != '#') {
            smiles_list.push_back(line);
        }
    }

    if (smiles_list.empty()) {
        throw std::runtime_error("No valid SMILES found in file: " + file_path);
    }

    return smiles_list;
}

std::vector<std::string> getMoleculeFiles(const std::string& dir_path) {
    std::vector<std::string> files;

    try {
        for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().string();
                std::string ext = filename.substr(filename.find_last_of('.') + 1);
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

                if (ext == "sdf") {
                    files.push_back(filename);
                }
            }
        }

        std::sort(files.begin(), files.end());

    } catch (const std::exception& e) {
        throw std::runtime_error("Cannot read directory " + dir_path + ": " + e.what());
    }

    return files;
}

Molecule loadMolecule(const std::string& source,
                      const SurfaceParams& surface_params) {
    if (std::filesystem::exists(source)) {
        std::string ext = source.substr(source.find_last_of('.') + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "sdf") {
            Molecule mol(source);
            return mol;
        } else if (ext == "txt") {
            // For text files, assume they contain SMILES strings and load the first one
            std::vector<std::string> smiles_list = loadSmilesFromFile(source);
            if (!smiles_list.empty()) {
                Molecule mol(smiles_list[0].c_str());
                if (mol.get_mol() != nullptr) {
                    
                    if (surface_params.conformer_params.enable_optimization) {
                        // Use advanced conformer generation
                        mol.generateAdvancedConformers(surface_params.conformer_params);
                    } else {
                        // Use simple conformer generation
                        mol.get_mol()->clearConformers();
                        mol.calculate_conformers(surface_params.conformer_params.num_clusters);
                    }
                }
                return mol;
            } else {
                throw std::runtime_error("No valid SMILES found in text file: " + source);
            }
        } else {
            throw std::runtime_error("Unsupported file format: " + ext +
                                     ". Only SDF and TXT files are supported.");
        }
    } else {
        // For SMILES strings, create molecule and regenerate conformers with specified number
        Molecule mol(source.c_str());
        if (mol.get_mol() != nullptr) {
            
            if (surface_params.conformer_params.enable_optimization) {
                // Use advanced conformer generation
                mol.generateAdvancedConformers(surface_params.conformer_params);
            } else {
                // Use simple conformer generation
                mol.calculate_conformers(surface_params.conformer_params.num_clusters);
            }
        }
        return mol;
    }
}

void computeProperty(Molecule& mol, int conformer_id, const std::string& property_name, const SurfaceParams& params) {
    try {
        // Convert CLI SurfaceParams to Molecule SurfaceParams
        SurfaceParams mol_params;
        mol_params.num_vertices = params.num_vertices;
        mol_params.radius = params.radius;
        mol_params.density = params.density;
        mol_params.hdensity = params.hdensity;
        mol_params.type = params.type;
        mol_params.program = params.program;
        mol_params.sample_method = params.sample_method;
        mol_params.properties = params.properties;
        mol_params.output_type = params.output_type;
        mol_params.include_normals = params.include_normals;
        mol_params.top_n = params.top_n;
        mol_params.conformer_params = params.conformer_params;
        mol_params.save_meshes = params.save_meshes;
        mol_params.charge_method = params.charge_method; // Pass the charge method
        
        mol.compute(property_name, conformer_id, mol_params);
        if (property_name == "esp") {
            // Get the ESP results
            std::vector<double> charges =
                mol.getProperty<std::vector<double>>("esp", "esp_charges", conformer_id);
            std::vector<double> potentials = mol.getProperty<std::vector<double>>(
                "esp", "electrostatic_potential", conformer_id);

            if (!potentials.empty()) {
#ifdef DEBUG
                double min_pot = *std::min_element(potentials.begin(), potentials.end());
                double max_pot = *std::max_element(potentials.begin(), potentials.end());
                double avg_pot =
                    std::accumulate(potentials.begin(), potentials.end(), 0.0) / potentials.size();
                std::cout << "ESP Potential Statistics:" << std::endl;
                std::cout << "  Actual min potential: " << std::fixed << std::setprecision(3)
                          << min_pot << " eV" << std::endl;
                std::cout << "  Actual max potential: " << std::fixed << std::setprecision(3)
                          << max_pot << " eV" << std::endl;
                std::cout << "  Max potential: " << std::fixed << std::setprecision(3) << max_pot
                          << " eV" << std::endl;
                std::cout << "  Average potential: " << std::fixed << std::setprecision(3)
                          << avg_pot << " eV" << std::endl;
                std::cout << "  Computed for " << potentials.size() << " surface points"
                          << std::endl;
#endif
            }

        } else if (property_name == "hydrophobicity") {
            // Get the hydrophobicity results
            auto hydrophobicity_values = mol.getProperty<std::vector<double>>(
                "hydrophobicity", "hydrophobicity_values", conformer_id);

            if (!hydrophobicity_values.empty()) {
#ifdef DEBUG
                auto min_hyd =
                    mol.getProperty<double>("hydrophobicity", "hydrophobicity_min", conformer_id);
                auto max_hyd =
                    mol.getProperty<double>("hydrophobicity", "hydrophobicity_max", conformer_id);
                auto avg_hyd =
                    mol.getProperty<double>("hydrophobicity", "hydrophobicity_avg", conformer_id);
                std::cout << "Hydrophobicity Statistics:" << std::endl;
                std::cout << "  Min hydrophobicity: " << std::fixed << std::setprecision(3)
                          << min_hyd << std::endl;
                std::cout << "  Max hydrophobicity: " << std::fixed << std::setprecision(3)
                          << max_hyd << std::endl;
                std::cout << "  Average hydrophobicity: " << std::fixed << std::setprecision(3)
                          << avg_hyd << std::endl;
                std::cout << "  Computed for " << hydrophobicity_values.size() << " surface points"
                          << std::endl;
#endif
            }

        } else if (property_name == "hbond") {
            // Get the hydrogen bond results
            auto hbond_values =
                mol.getProperty<std::vector<double>>("hbond", "hbond_values", conformer_id);

            if (!hbond_values.empty()) {
#ifdef DEBUG
                auto min_hbond = mol.getProperty<double>("hbond", "hbond_min", conformer_id);
                auto max_hbond = mol.getProperty<double>("hbond", "hbond_max", conformer_id);
                auto avg_hbond = mol.getProperty<double>("hbond", "hbond_avg", conformer_id);
                std::cout << "Hydrogen Bond Statistics:" << std::endl;
                std::cout << "  Min hbond potential: " << std::fixed << std::setprecision(3)
                          << min_hbond << std::endl;
                std::cout << "  Max hbond potential: " << std::fixed << std::setprecision(3)
                          << max_hbond << std::endl;
                std::cout << "  Average hbond potential: " << std::fixed << std::setprecision(3)
                          << avg_hbond << std::endl;
                std::cout << "  Computed for " << hbond_values.size() << " surface points"
                          << std::endl;
#endif
            }
        }

    } catch (const std::exception& e) {
        std::cout << "  ✗ " << property_name << " calculation failed: " << e.what() << std::endl;
    }
}

void alignMolecule(Molecule& reference_mol, Molecule& input_mol, const std::string& output_dir,
                   const SurfaceParams& params, const std::string& input_name) {
    // Create reference surface
    int ref_conformer_id =
        reference_mol.createSurface(0, params.num_vertices, params.radius, params.density,
                                    params.hdensity, params.type, params.program, params.sample_method, params.addH);

    // Check if this is a SMILES molecule (multiple conformers)
    bool is_smiles_molecule = input_mol.getNumConformers() > 1;

    if ((is_smiles_molecule && params.conformer_params.num_clusters > 1) || params.conformer_params.enable_optimization) {
        // Multi-conformer alignment for SMILES molecules
        std::cout << "\n=== Multi-Conformer Alignment ===" << std::endl;
        std::cout << "Processing " << input_mol.getNumConformers() << " conformers..." << std::endl;

        std::vector<std::pair<int, Molecule::AlignmentResult>> conformer_results;

        // Process each conformer
        for (int conf_id = 0; conf_id < static_cast<int>(input_mol.getNumConformers()); ++conf_id) {

            // Create surface for this conformer
            int input_conformer_id =
                input_mol.createSurface(conf_id, params.num_vertices, params.radius, params.density,
                                        params.hdensity, params.type, params.program, params.sample_method, params.addH);

            // Compute requested surface properties BEFORE alignment (for both molecules)
            if (!params.properties.empty()) {
                for (const auto& property : params.properties) {
                    computeProperty(reference_mol, ref_conformer_id, property, params);
                }

                for (const auto& property : params.properties) {
                    computeProperty(input_mol, input_conformer_id, property, params);
                }
            }

            // Perform alignment and get detailed results
            auto alignment_result = input_mol.align_to_molecule_detailed(
                reference_mol, ref_conformer_id, input_conformer_id);

            // Apply transformation to surface FIRST (before conformer transformation invalidates
            // cache)
            input_mol.transformSurface(input_conformer_id, alignment_result.transformation);

            // Save individual mesh files with temporary names (BEFORE conformer transformation invalidates surface cache)
            if (params.save_meshes) {
                std::string temp_surface_filename = output_dir + "/" + input_name + "_temp_conf" +
                                                   std::to_string(conf_id) + "_surface.ply";
                if (params.output_type == "mesh") {
                    input_mol.saveSurfaceMesh(input_conformer_id, temp_surface_filename,
                                              params.include_normals, true, params.properties);
                } else {
                    input_mol.saveSurfacePointCloud(input_conformer_id, temp_surface_filename,
                                                    params.include_normals, true, params.properties);
                }
            }

            // Apply transformation to conformer coordinates (this invalidates surface cache)
            input_mol.transformConformer(input_conformer_id, alignment_result.transformation);

            // Store results
            conformer_results.push_back({conf_id, alignment_result});

            // std::cout << "  Conformer " << conf_id << " alignment:" << std::endl;
            // std::cout << "    Fitness: " << std::fixed << std::setprecision(6)
            //           << alignment_result.fitness << std::endl;
            // std::cout << "    RMSE: " << std::fixed << std::setprecision(6)
            //           << alignment_result.inlier_rmse << std::endl;
            // std::cout << "    Used Colored ICP: "
                    //   << (alignment_result.used_colored_icp ? "Yes" : "No") << std::endl;
        }

        // Sort conformers by fitness (higher is better) - ALWAYS sort regardless of top_n
        std::sort(conformer_results.begin(), conformer_results.end(),
                 [](const std::pair<int, Molecule::AlignmentResult>& a,
                    const std::pair<int, Molecule::AlignmentResult>& b) {
                     return a.second.fitness > b.second.fitness;  // Higher fitness is better
                 });
        
        // Keep only top N conformers if specified
        if (params.top_n > 0) {
            if (static_cast<int>(conformer_results.size()) > params.top_n) {
                conformer_results.resize(params.top_n);
                std::cout << "\n✓ Kept top " << params.top_n << " conformers out of " 
                          << input_mol.getNumConformers() << " total" << std::endl;
            }
        }

        // Rename individual surface files with rank and fitness in filename (if requested)
        if (params.save_meshes) {
            for (size_t rank = 0; rank < conformer_results.size(); ++rank) {
                int conf_id = conformer_results[rank].first;
                double fitness = conformer_results[rank].second.fitness;
                
                // Create new filename with rank and fitness: e.g., "imatinib_simi_rank1_f841.ply"
                std::string new_surface_filename = output_dir + "/" + input_name + "_rank" + 
                                                 std::to_string(rank + 1) + "_f" + 
                                                 std::to_string(static_cast<int>(fitness * 1000)) + ".ply";
                
                // Rename the temporary file
                std::string temp_surface_filename = output_dir + "/" + input_name + "_temp_conf" +
                                                   std::to_string(conf_id) + "_surface.ply";
                
                if (std::rename(temp_surface_filename.c_str(), new_surface_filename.c_str()) != 0) {
                    std::cerr << "Warning: Could not rename " << temp_surface_filename << " to " << new_surface_filename << std::endl;
                }
            }
            
            // Clean up any remaining temporary files for conformers not in top N
            for (int conf_id = 0; conf_id < static_cast<int>(input_mol.getNumConformers()); ++conf_id) {
                std::string temp_surface_filename = output_dir + "/" + input_name + "_temp_conf" +
                                                   std::to_string(conf_id) + "_surface.ply";
                std::remove(temp_surface_filename.c_str());
            }
        }

        // Save all conformers with alignment scores to a single SDF file
        std::string conformer_filename =
            output_dir + "/" + input_name + "_aligned.sdf";
        input_mol.saveMultiConformerSDF(conformer_results, conformer_filename, true);
        std::cout << "\n✓ Multi-conformer SDF saved: " << conformer_filename << std::endl;

        // Save reference surface (only once)
        if (params.output_type == "mesh") {
            reference_mol.saveSurfaceMesh(ref_conformer_id, output_dir + "/reference_surface.ply",
                                          params.include_normals, true, params.properties);
        } else {
            reference_mol.saveSurfacePointCloud(ref_conformer_id,
                                                output_dir + "/reference_surface.ply",
                                                params.include_normals, true, params.properties);
        }

    } else {
        // Single conformer alignment (original behavior)
        int input_conformer_id =
            input_mol.createSurface(0, params.num_vertices, params.radius, params.density,
                                    params.hdensity, params.type, params.program, params.sample_method, params.addH);

        // Compute requested surface properties BEFORE alignment (for both molecules)
        if (!params.properties.empty()) {
            for (const auto& property : params.properties) {
                computeProperty(reference_mol, ref_conformer_id, property, params);
            }

            for (const auto& property : params.properties) {
                computeProperty(input_mol, input_conformer_id, property, params);
            }
        } else {
            std::cout << "\nSkipping surface property computation (--properties none)" << std::endl;
        }

// Perform alignment (uses already computed properties for coloring)
#ifdef DEBUG
        std::cout << "Aligning " << input_mol.get_smiles() << " to " << reference_mol.get_smiles()
                  << std::endl;
#endif
        Eigen::Isometry3d transformation =
            input_mol.align_to_molecule(reference_mol, ref_conformer_id, input_conformer_id);

        // Apply transformation to surface FIRST (before conformer transformation invalidates cache)
        input_mol.transformSurface(input_conformer_id, transformation);

        // Save aligned conformer and surface
        std::string conformer_filename = output_dir + "/" + input_name + "_aligned.sdf";
        std::string surface_filename = output_dir + "/" + input_name + "_surface.ply";

        // Save surface based on output type (BEFORE conformer transformation invalidates surface
        // cache) - only if --save-meshes is specified
        if (params.save_meshes) {
            if (params.output_type == "mesh") {
                input_mol.saveSurfaceMesh(input_conformer_id, surface_filename, params.include_normals,
                                          true, params.properties);
                reference_mol.saveSurfaceMesh(ref_conformer_id, output_dir + "/reference_surface.ply",
                                              params.include_normals, true, params.properties);
            } else {
                input_mol.saveSurfacePointCloud(input_conformer_id, surface_filename,
                                                params.include_normals, true, params.properties);
                reference_mol.saveSurfacePointCloud(ref_conformer_id,
                                                    output_dir + "/reference_surface.ply",
                                                    params.include_normals, true, params.properties);
            }
        }

        // Apply transformation to conformer coordinates (this invalidates surface cache)
        input_mol.transformConformer(input_conformer_id, transformation);

        // Save conformer SDF (after conformer coordinates are transformed)
        input_mol.saveConformerSDF(input_conformer_id, conformer_filename, true);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printUsage(argv[0]);
        return 1;
    }

    std::string reference_source = argv[1];
    std::string input_source = argv[2];
    std::string output_dir = argv[3];

    int arg_index = 4;
    SurfaceParams surface_params = parseSurfaceParams(argc, argv, arg_index);

    try {
        std::filesystem::create_directories(output_dir);
#ifdef DEBUG
        std::cout << "Loading reference molecule from: " << reference_source << std::endl;
#endif
        Molecule reference_mol = loadMolecule(reference_source, surface_params);
#ifdef DEBUG
        std::cout << "Reference molecule: " << reference_mol.get_smiles() << std::endl;
#endif
        reference_mol.saveConformerSDF(0, output_dir + "/reference_conformer.sdf", true);
        reference_mol.createSurface(0, surface_params.num_vertices, surface_params.radius,
                                    surface_params.density, surface_params.hdensity,
                                    surface_params.type, surface_params.program, surface_params.sample_method, surface_params.addH);

        if (std::filesystem::is_directory(input_source)) {
#ifdef DEBUG
            std::cout << "\nProcessing directory: " << input_source << std::endl;
#endif

            std::vector<std::string> files = getMoleculeFiles(input_source);
#ifdef DEBUG
            std::cout << "Found " << files.size() << " SDF file(s)" << std::endl;
#endif

            for (size_t i = 0; i < files.size(); ++i) {
                const std::string& file_path = files[i];
                std::string filename = std::filesystem::path(file_path).stem().string();
#ifdef DEBUG
                std::cout << "\nProcessing file " << (i + 1) << "/" << files.size() << ": "
                          << filename << std::endl;
#endif
                try {
                    // Load molecules from file
                    auto molecules = Molecule::loadFromMultiMolFile(file_path);
#ifdef DEBUG
                    std::cout << "  Loaded " << molecules.size() << " molecule(s) from file"
                              << std::endl;
#endif

                    // Align each molecule
                    for (size_t j = 0; j < molecules.size(); ++j) {
                        std::string mol_name = filename + "_mol_" + std::to_string(j + 1);
                        alignMolecule(reference_mol, *molecules[j], output_dir, surface_params,
                                      mol_name);
                    }

                } catch (const std::exception& e) {
                    std::cerr << "  ✗ Error processing file " << filename << ": " << e.what()
                              << std::endl;
                }
            }

        } else if (std::filesystem::exists(input_source)) {
            // Check if it's a text file with SMILES
            std::string ext = input_source.substr(input_source.find_last_of('.') + 1);
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

            if (ext == "txt") {
// Process text file with multiple SMILES
#ifdef DEBUG
                std::cout << "\nProcessing SMILES text file: " << input_source << std::endl;
#endif

                try {
                    std::vector<std::string> smiles_list = loadSmilesFromFile(input_source);
                    std::string base_name = std::filesystem::path(input_source).stem().string();

#ifdef DEBUG
                    std::cout << "Found " << smiles_list.size() << " SMILES string(s)" << std::endl;
#endif

                    for (size_t i = 0; i < smiles_list.size(); ++i) {
                        std::cout << "\nProcessing SMILES " << (i + 1) << "/" << smiles_list.size()
                                  << ": " << smiles_list[i] << std::endl;

                        try {
                            // Use loadMolecule to properly handle both simple and advanced
                            // conformer generation
                            Molecule input_mol =
                                loadMolecule(smiles_list[i], surface_params);

                            std::string mol_name = base_name + "_smiles_" + std::to_string(i + 1);
                            alignMolecule(reference_mol, input_mol, output_dir, surface_params,
                                          mol_name);

                        } catch (const std::exception& e) {
                            std::cerr << "  ✗ Error processing SMILES " << (i + 1) << ": "
                                      << e.what() << std::endl;
                        }
                    }

                } catch (const std::exception& e) {
                    std::cerr << "✗ Error processing text file: " << e.what() << std::endl;
                    return 1;
                }
            } else {
                // Process as SDF file (original logic)
                try {
                    // For single files, try single molecule loading first
                    try {
                        Molecule input_mol(input_source);
                        std::string base_name = std::filesystem::path(input_source).stem().string();
                        alignMolecule(reference_mol, input_mol, output_dir, surface_params,
                                      base_name);
                    } catch (const std::exception& e) {
                        // Fall back to multi-molecule loader
                        auto molecules = Molecule::loadFromMultiMolFile(input_source);

                        // Align each molecule
                        std::string base_name = std::filesystem::path(input_source).stem().string();
                        for (size_t i = 0; i < molecules.size(); ++i) {
                            std::string mol_name = base_name + "_mol_" + std::to_string(i + 1);
                            alignMolecule(reference_mol, *molecules[i], output_dir, surface_params,
                                          mol_name);
                        }
                    }

                } catch (const std::exception& e) {
                    std::cerr << "✗ Error processing file: " << e.what() << std::endl;
                    return 1;
                }
            }

        } else {
// Assume it's a SMILES string
#ifdef DEBUG
            std::cout << "\nProcessing SMILES: " << input_source << std::endl;
#endif
            try {
                Molecule input_mol = loadMolecule(input_source, surface_params);
                alignMolecule(reference_mol, input_mol, output_dir, surface_params, "smiles");

            } catch (const std::exception& e) {
                std::cerr << "✗ Error processing SMILES: " << e.what() << std::endl;
                return 1;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "✓ Alignment completed successfully!" << std::endl;
    std::cout << "Results saved to: " << output_dir << std::endl;
    return 0;
}