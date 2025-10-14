#include "surface.h"
#include "fps.h"
#include "props/esp.h"
#include "props/property_params.h"
#include "utils.h"
#include "surface_fibonacci.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <random>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <open3d/Open3D.h>

using PointCloudData = std::vector<Eigen::Vector3f>;

// Van der Waals radii (Grimme) in Angstroms
const std::map<std::string, double> Surface::vdw_radii_ = {
    {"H", 1.09155},  {"He", 0.86735}, {"Li", 1.7478},  {"Be", 1.5491},  {"B", 1.608},
    {"C", 1.45515},  {"N", 1.31125},  {"O", 1.24085},  {"F", 1.1498},   {"Ne", 1.0687},
    {"Na", 1.8541},  {"Mg", 1.74195}, {"Al", 2.0053},  {"Si", 1.89585}, {"P", 1.75085},
    {"S", 1.65535},  {"Cl", 1.5523},  {"Ar", 1.4574},  {"K", 2.12055},  {"Ca", 2.05175},
    {"Sc", 1.94515}, {"Ti", 1.8821},  {"V", 1.86055},  {"Cr", 1.7207},  {"Mn", 1.7731},
    {"Fe", 1.72105}, {"Co", 1.71635}, {"Ni", 1.6731},  {"Cu", 1.6504},  {"Zn", 1.61545},
    {"Ga", 1.97895}, {"Ge", 1.93095}, {"As", 1.83125}, {"Se", 1.7634},  {"Br", 1.6831},
    {"Kr", 1.6048},  {"Rb", 2.3088},  {"Sr", 2.2382},  {"Y", 2.1098},   {"Zr", 2.02985},
    {"Nb", 1.9298},  {"Mo", 1.87715}, {"Tc", 1.7845},  {"Ru", 1.73115}, {"Rh", 1.69875},
    {"Pd", 1.67625}, {"Ag", 1.6654},  {"Cd", 1.731},   {"In", 2.13115}, {"Sn", 2.0937},
    {"Sb", 2.0075},  {"Te", 1.94505}, {"I", 1.869},    {"Xe", 1.79445}, {"Cs", 2.52835},
    {"Ba", 2.5907},  {"La", 2.31305}, {"Ce", 2.31005}, {"Pr", 2.2851},  {"Nd", 2.26355},
    {"Pm", 2.2448},  {"Sm", 2.22575}, {"Eu", 2.2117},  {"Gd", 2.06215}, {"Tb", 2.12135},
    {"Dy", 2.07705}, {"Ho", 2.1397},  {"Er", 2.1225},  {"Tm", 2.1104},  {"Yb", 2.0993},
    {"Lu", 2.0065},  {"Hf", 2.1225},  {"Ta", 2.049},   {"W", 1.99275},  {"Re", 1.94775},
    {"Os", 1.8745},  {"Ir", 1.7228},  {"Pt", 1.67625}, {"Au", 1.6282},  {"Hg", 1.67995},
    {"Tl", 2.15635}, {"Pb", 2.1382},  {"Bi", 2.05875}, {"Po", 2.0027},  {"At", 1.9322},
    {"Rn", 1.8608},  {"Fr", 2.5398},  {"Ra", 2.4647},  {"Ac", 2.35215}, {"Th", 2.2126},
    {"Pa", 2.2297},  {"U", 2.19785},  {"Np", 2.17695}, {"Pu", 2.21705}};

Surface::Surface(const Molecule& molecule, int conformer_id, int num_vertices, double radius,
                 double density, double hdensity, std::string type, std::string program,
                 std::string sample_method)
    : molecule(molecule), conformer_id(conformer_id) {
    assert(num_vertices > 4 && "num_vertices must be greater than 4");
    assert(radius > 0 && "radius must be greater than 0");
    assert(density > 0 && "density must be greater than 0");
    assert(hdensity > 0 && "hdensity must be greater than 0");
    assert((type == "tses" || type == "ases") && "type must be either tses or ases");
    assert((program == "msms" || program == "fibonacci") && "program must be either msms or fibonacci");
    assert((sample_method == "random" || sample_method == "fps" || sample_method == "full") &&
           "sample_method must be random, fps, or full");

    // Validate conformer ID using public method
    if (molecule.getNumConformers() <= conformer_id) {
        throw std::runtime_error("Conformer " + std::to_string(conformer_id) + " does not exist");
    }

    this->num_vertices = num_vertices;
    this->radius = radius;
    this->density = density;
    this->hdensity = hdensity;
    this->type = type;
    this->sample_method = sample_method;

    // Use requested program; if unavailable, fall back to fibonacci (no external deps)
    if (program == "fibonacci") {
        this->program = program;
    } else if (program == "msms" && !isAvailable("msms")) {
#ifdef DEBUG
        std::cerr << "MSMS is not available; falling back to fibonacci" << std::endl;
#endif
        this->program = "fibonacci";
    } else {
        this->program = program;
    }

    // Automatically generate surface upon construction
    generateSurface("molecule", "surface");

    // Apply subsampling based on sample_method
    applySubsampling();
}

Surface::~Surface() {
}

void Surface::generateSurface(std::string_view input_file, std::string_view output_file) {
    // Use Fibonacci if requested or if MSMS is unavailable
    if (program == "fibonacci" || !isAvailable("msms")) {
        try {
            generateFibonacciSurface();
        } catch (const std::exception& e) {
            std::cerr << "Error generating Fibonacci surface: " << e.what() << std::endl;
            throw;
        } catch (...) {
            std::cerr << "Error generating Fibonacci surface: Unknown error" << std::endl;
            throw;
        }
        return;
    }
    
    std::string input_str(input_file);
    std::string extension = getFileExtension(input_file);

    std::string xyzr_file;

    if (extension == "xyzr") {
        xyzr_file = input_str;
        generateSurfaceFromXYZR(xyzr_file, output_file);
    } else {
        xyzr_file = std::string(output_file) + "_temp.xyzr";
        generate3DCoordinates();

        generateXYZRFile(xyzr_file);

        generateSurfaceFromXYZR(xyzr_file, output_file);
        safe_remove(xyzr_file);
    }
}

void Surface::generateSurfaceFromXYZR(std::string_view xyzr_file, std::string_view output_file) {
    std::ifstream in{std::string(xyzr_file)};
    if (!in) {
        throw std::runtime_error("Failed to open xyzr file: " + std::string(xyzr_file));
    }
    int line_count = 0;
    std::string line;
    while (std::getline(in, line)) {
        line_count++;
        auto columns = split_on_whitespace(line);
        if (columns.size() < EXPECTED_MIN_COLUMNS) {
            throw std::runtime_error("Line " + std::to_string(line_count) + " has less than " +
                                     std::to_string(EXPECTED_MIN_COLUMNS) + " columns");
        }
    }

    if (program == "msms") {
        run_msms(xyzr_file, output_file);
    } else {
        throw std::runtime_error("Unsupported program: " + program);
    }
}

void Surface::run_msms(std::string_view xyzr_file, std::string_view output_file) {
    const std::string if_arg = std::string(xyzr_file);
    const std::string of_arg = std::string(output_file);
    std::string command = program + " -probe_radius " + std::to_string(radius) + " -density " +
                          std::to_string(density) + " -hdensity " + std::to_string(hdensity) +
                          " -if " + if_arg + " -of " + of_arg;
#ifndef DEBUG
    command += " > /dev/null 2>&1";
#endif
    if (!run_command(command)) {
        throw std::runtime_error("Failed to run " + program + " command: " + command);
    }
    read_vert_file(of_arg + ".vert");
    read_face_file(of_arg + ".face");

    safe_remove(of_arg + ".vert");
    safe_remove(of_arg + ".face");
}

bool Surface::read_vert_file(std::string_view vert_file) {
    std::ifstream in{std::string(vert_file)};
    if (!in) {
        throw std::runtime_error("Failed to open vert file: " + std::string(vert_file));
    }

    std::string line;
    int vertex_index = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        auto columns = split_on_whitespace(line);
        if (columns.size() == 0) {
            continue;
        }
        if (columns.size() == 4) {
            continue;
        }
        if (columns.size() < EXPECTED_VERT_COLUMNS) {
            throw std::runtime_error(
                "Line in .vert file has wrong number of columns. Expected >= " +
                std::to_string(EXPECTED_VERT_COLUMNS) + ", got " + std::to_string(columns.size()));
        }
        const std::string& key = columns.at(9);
        this->vertices[vertex_index] = {key,
                                        {
                                            std::stod(columns.at(0)),
                                            std::stod(columns.at(1)),
                                            std::stod(columns.at(2)),
                                            std::stod(columns.at(3)),
                                            std::stod(columns.at(4)),
                                            std::stod(columns.at(5)),
                                        }};
        vertex_index++;
    }
    return true;
}

void Surface::applySubsampling() {
    if (sample_method == "full") {
        return;
    }

    size_t current_vertices = vertices.size();
    if (current_vertices <= static_cast<size_t>(num_vertices)) {
        // Already have fewer or equal vertices than requested
        return;
    }
    
    // Warn user that vertices will be subsampled
    std::cout << "Warning: Surface has " << current_vertices << " vertices but --vertices=" 
              << num_vertices << " requested. Subsampling using " << sample_method << " method." << std::endl;

    std::map<int, std::pair<std::string, std::vector<double>>> sampled_vertices;

    if (sample_method == "random") {
        sampled_vertices = select_n_vertices_random(num_vertices);
    } else if (sample_method == "fps") {
        sampled_vertices = select_n_vertices_fps(num_vertices);
    } else {
        throw std::runtime_error("Unknown sampling method: " + sample_method);
    }

    // Replace the original vertices with sampled ones
    vertices = sampled_vertices;
}

bool Surface::read_face_file(std::string_view face_file) {
    std::ifstream in{std::string(face_file)};
    if (!in) {
        throw std::runtime_error("Failed to open face file: " + std::string(face_file));
    }

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            continue;  // skip comment lines
        }
        auto columns = split_on_whitespace(line);
        if (columns.empty()) {
            continue;
        }

        if (columns.size() != 5) {
            continue;
        }
        // only the first three columns are the vertices connected by the face
        // we need to subtract 1 because the vertices are 0-indexed
        this->faces.push_back({std::stoi(columns.at(0)) - 1, std::stoi(columns.at(1)) - 1,
                               std::stoi(columns.at(2)) - 1});
    }
    return true;
}

int Surface::get_num_vertices() const {
    return vertices.size();
}

std::map<int, std::pair<std::string, std::vector<double>>> Surface::select_n_vertices(
    size_t n_vertices, std::string sampling_method) {
    if (n_vertices > vertices.size()) {
        throw std::runtime_error("n_vertices is greater than the number of vertices in the file");
    }
    if (sampling_method == "random") {
        return select_n_vertices_random(n_vertices);
    } else if (sampling_method == "fps") {
        return select_n_vertices_fps(n_vertices);
    } else {
        throw std::runtime_error("Invalid sampling method: " + sampling_method);
    }
}

std::map<int, std::pair<std::string, std::vector<double>>> Surface::select_n_vertices_random(
    size_t n_vertices) {
    static const int numVertices = vertices.size();
    std::vector<int> v(numVertices + 1);
    for (size_t i = 0; i < v.size(); ++i)
        v[i] = i;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(v.begin(), v.end(), gen);
    std::map<int, std::pair<std::string, std::vector<double>>> selected_vertices;

    for (size_t i = 0; i < n_vertices; i++) {
        // Safety check: ensure index is within bounds
        if (v[i] < 0 || static_cast<size_t>(v[i]) >= vertices.size()) {
            std::cerr << "Warning: Invalid vertex index " << v[i] 
                      << " (total vertices: " << vertices.size() 
                      << "), skipping vertex " << i << std::endl;
            continue;
        }
        selected_vertices[i] = vertices[v[i]];
    }
    return selected_vertices;
}

std::map<int, std::pair<std::string, std::vector<double>>> Surface::select_n_vertices_fps(
    size_t n_vertices) {
    if (n_vertices > vertices.size()) {
        throw std::runtime_error("n_vertices is greater than the number of vertices in the file");
    }
    std::map<int, std::pair<std::string, std::vector<double>>> selected_vertices;
    std::vector<Point> points;
    for (const auto& [index, vertex] : vertices) {
        points.push_back({vertex.second[0], vertex.second[1], vertex.second[2]});
    }
    std::vector<int> sampled_indices = farthestPointSampling(points, n_vertices);
    for (size_t i = 0; i < sampled_indices.size(); i++) {
        // Safety check: ensure index is within bounds
        if (sampled_indices[i] < 0 || static_cast<size_t>(sampled_indices[i]) >= vertices.size()) {
            std::cerr << "Warning: Invalid sampled vertex index " << sampled_indices[i] 
                      << " (total vertices: " << vertices.size() 
                      << "), skipping vertex " << i << std::endl;
            continue;
        }
        
        const auto& vertex = vertices[sampled_indices[i]];
        // Safety check: ensure vertex has enough coordinates
        if (vertex.second.size() < 6) {
            std::cerr << "Warning: Vertex " << sampled_indices[i] 
                      << " has insufficient coordinates (" << vertex.second.size() 
                      << "), skipping" << std::endl;
            continue;
        }
        
        selected_vertices[i] = {
            vertex.first,
            {vertex.second[0], vertex.second[1], vertex.second[2],
             vertex.second[3], vertex.second[4], vertex.second[5]}};
    }
    return selected_vertices;
}

void Surface::registerProperty(std::unique_ptr<ISurfaceProperty> prop) {
    props_[prop->key()] = std::move(prop);
}

void Surface::computeProperty(const std::string& key) {
    if (cache_.count(key))
        return;
    auto it = props_.find(key);
    if (it == props_.end())
        throw std::runtime_error("Unknown property: " + key);
    it->second->compute(*this, cache_);
}

template <typename T>
const T& Surface::getProperty(const std::string& key) const {
    auto it = cache_.find(key);
    if (it == cache_.end())
        throw std::runtime_error("Property not computed: " + key);
    return std::any_cast<const T&>(it->second);
}

void Surface::invalidateProperties() {
    cache_.clear();
}

double Surface::getVertexColorValue(
    const std::string& property_name, int vertex_index, const std::string& fallback_atom_name,
    const std::unordered_map<std::string, std::any>* property_cache) const {
    (void)fallback_atom_name; // Suppress unused parameter warning
    // If no property specified, return neutral value (0.0)
    if (property_name.empty()) {
        return 0.0;
    }

    // Try to get property values from cache
    if (property_cache != nullptr) {
        try {
            std::string cache_key;

            // Dynamic property name to cache key mapping
            std::vector<std::string> possible_keys = {
                property_name + "_values",     // e.g., "esp_values"
                "electrostatic_potential",     // ESP specific
                property_name + "_potential",  // e.g., "hbond_potential"
                property_name + "_charges",    // e.g., "esp_charges"
                property_name                  // Direct property name
            };

            // Find the first matching key in the cache
            for (const auto& key : possible_keys) {
                if (property_cache->find(key) != property_cache->end()) {
                    cache_key = key;
                    break;
                }
            }

            if (cache_key.empty()) {
                cache_key = property_name + "_values";  // Fallback
            }

            // Check if property exists in cache
            auto it = property_cache->find(cache_key);
            if (it != property_cache->end()) {
                // Get the property values vector
                const auto& values = std::any_cast<const std::vector<double>&>(it->second);

                if (vertex_index < static_cast<int>(values.size())) {
                    double value = values[vertex_index];
                    return normalizePropertyValue(property_name, value, property_cache);
                }
            }
        } catch (const std::exception& e) {
            // Property not found or wrong type, fall back to neutral value
        }
    }

    // Fallback to neutral value
    return 0.0;
}

double Surface::normalizePropertyValue(
    const std::string& property_name, double value,
    const std::unordered_map<std::string, std::any>* property_cache) const {
    // Normalize value to [-1, 1] range using property-specific min/max
    try {
        // First try to get min/max from PropertyParamsRegistry
        PropertyParams params = PropertyParamsRegistry::getParams(property_name);
        double min_val = params.min_value;
        double max_val = params.max_value;

        if (max_val > min_val) {
            // Normalize to [0, 1] first
            double normalized = (value - min_val) / (max_val - min_val);
            normalized = std::max(0.0, std::min(1.0, normalized));  // Clamp to [0, 1]

            // Convert to [-1, 1] range
            return 2.0 * normalized - 1.0;
        }
    } catch (const std::exception&) {
        // Property not found in registry, try cache-based approach
        try {
            std::string min_key = property_name + "_min";
            std::string max_key = property_name + "_max";

            // Get min and max values for this property
            if (property_cache != nullptr) {
                auto min_it = property_cache->find(min_key);
                auto max_it = property_cache->find(max_key);

                if (min_it != property_cache->end() && max_it != property_cache->end()) {
                    double min_val = std::any_cast<double>(min_it->second);
                    double max_val = std::any_cast<double>(max_it->second);

                    if (max_val > min_val) {
                        // Normalize to [0, 1] first
                        double normalized = (value - min_val) / (max_val - min_val);
                        normalized = std::max(0.0, std::min(1.0, normalized));  // Clamp to [0, 1]

                        // Convert to [-1, 1] range
                        return 2.0 * normalized - 1.0;
                    }
                }
            }
        } catch (const std::exception&) {
            // Use default range as last resort
            double min_val = -1.0, max_val = 1.0;
            double normalized = (value - min_val) / (max_val - min_val);
            normalized = std::max(0.0, std::min(1.0, normalized));
            return 2.0 * normalized - 1.0;
        }
    }

    return 0.0;  // Neutral value
}

std::tuple<int, int, int> Surface::getVertexColor(
    const std::string& property_name, int vertex_index, const std::string& fallback_atom_name,
    const std::unordered_map<std::string, std::any>* property_cache) const {
    // If no property specified, use atom type coloring
    if (property_name.empty()) {
        return getAtomTypeColor(fallback_atom_name);
    }

    // Try to get property values from cache
    if (property_cache != nullptr) {
        try {
            std::string cache_key;

            // Map property names to their cache keys
            if (property_name == "esp") {
                cache_key = "electrostatic_potential";
            } else if (property_name == "hydrophobicity") {
                cache_key = "hydrophobicity_values";
            } else if (property_name == "hbond") {
                cache_key = "hbond_values";
            } else {
                // Try generic "_values" suffix
                cache_key = property_name + "_values";
            }

            // Check if property exists in cache
            auto it = property_cache->find(cache_key);
            if (it != property_cache->end()) {
                // Get the property values vector
                const auto& values = std::any_cast<const std::vector<double>&>(it->second);

                if (vertex_index < static_cast<int>(values.size())) {
                    double value = values[vertex_index];
                    return propertyValueToColor(property_name, value, property_cache);
                }
            }
        } catch (const std::exception& e) {
            // Property not found or wrong type, fall back to atom type coloring
        }
    }

    // Fallback to atom type coloring
    return getAtomTypeColor(fallback_atom_name);
}

std::tuple<int, int, int> Surface::getAtomTypeColor(const std::string& atom_name) const {
    // Standard CPK coloring scheme
    if (atom_name == "C")
        return {144, 144, 144};  // Carbon: dark gray
    else if (atom_name == "O")
        return {255, 0, 0};  // Oxygen: red
    else if (atom_name == "N")
        return {0, 0, 255};  // Nitrogen: blue
    else if (atom_name == "H")
        return {255, 255, 255};  // Hydrogen: white
    else if (atom_name == "S")
        return {255, 255, 0};  // Sulfur: yellow
    else if (atom_name == "P")
        return {255, 128, 0};  // Phosphorus: orange
    else if (atom_name == "F")
        return {0, 255, 0};  // Fluorine: green
    else if (atom_name == "Cl")
        return {0, 255, 0};  // Chlorine: green
    else if (atom_name == "Br")
        return {165, 42, 42};  // Bromine: brown
    else if (atom_name == "I")
        return {148, 0, 211};  // Iodine: purple
    else
        return {128, 128, 128};  // Default: gray
}

std::tuple<int, int, int> Surface::propertyValueToColor(
    const std::string& property_name, double value,
    const std::unordered_map<std::string, std::any>* property_cache) const {
    // Normalize value to [0, 1] range using property-specific min/max
    double normalized_value = 0.5;  // Default middle value

    try {
        std::string min_key, max_key;

        // Map property names to their cache keys
        if (property_name == "esp") {
            min_key = "electrostatic_potential_min";
            max_key = "electrostatic_potential_max";
        } else if (property_name == "hydrophobicity") {
            min_key = "hydrophobicity_min";
            max_key = "hydrophobicity_max";
        } else if (property_name == "hbond") {
            min_key = "hbond_min";
            max_key = "hbond_max";
        } else {
            // Try generic suffixes
            min_key = property_name + "_min";
            max_key = property_name + "_max";
        }

        // Get min and max values for this property
        if (property_cache != nullptr) {
            auto min_it = property_cache->find(min_key);
            auto max_it = property_cache->find(max_key);

            if (min_it != property_cache->end() && max_it != property_cache->end()) {
                double min_val = std::any_cast<double>(min_it->second);
                double max_val = std::any_cast<double>(max_it->second);

                if (max_val > min_val) {
                    normalized_value = (value - min_val) / (max_val - min_val);
                    normalized_value =
                        std::max(0.0, std::min(1.0, normalized_value));  // Clamp to [0, 1]
                }
            }
        }
    } catch (const std::exception& e) {
        // Use default middle value if min/max not available
    }

    // Universal color mapping: Min value → Red, Max value → Blue
    // This ensures consistent color interpretation across all properties
    int red = static_cast<int>(255 * (1.0 - normalized_value));  // Red for min values
    int blue = static_cast<int>(255 * normalized_value);         // Blue for max values

    // Clamp values to valid RGB range
    red = std::max(0, std::min(255, red));
    blue = std::max(0, std::min(255, blue));

    return {red, 0, blue};
}

std::string Surface::getFileExtension(std::string_view filename) const {
    std::string filename_str(filename);
    size_t dot_pos = filename_str.find_last_of('.');
    if (dot_pos == std::string::npos) {
        return "";  // No extension
    }
    return filename_str.substr(dot_pos + 1);
}

void Surface::generate3DCoordinates() {
    auto rdkit_mol = molecule.get_mol();
    if (!rdkit_mol) {
        throw std::runtime_error("Failed to get RDKit molecule for coordinate generation");
    }

    bool has_3d = false;
    if (rdkit_mol->getNumConformers() > 0) {
        auto conf = rdkit_mol->getConformer(this->conformer_id);
        if (conf.is3D()) {
            has_3d = true;
        }
    }

    if (!has_3d) {
        const_cast<Molecule&>(molecule).calculate_conformers(1);
    }
}

void Surface::generateXYZRFile(std::string_view output_file) {
    // Molecule reference is always valid (checked at construction)

    auto rdkit_mol = molecule.get_mol();
    if (!rdkit_mol) {
        throw std::runtime_error("Failed to get RDKit molecule for XYZR generation");
    }

    if (rdkit_mol->getNumConformers() == 0) {
        throw std::runtime_error("Molecule has no conformers - generate 3D coordinates first");
    }

    auto conf = rdkit_mol->getConformer(this->conformer_id);

    std::ofstream out{std::string(output_file)};
    if (!out) {
        throw std::runtime_error("Failed to create XYZR file: " + std::string(output_file));
    }

    for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
        auto atom = rdkit_mol->getAtomWithIdx(i);
        auto pos = conf.getAtomPos(i);

        std::string symbol = atom->getSymbol();
        double vdw_radius = DEFAULT_VDW_RADIUS;
        auto it = vdw_radii_.find(symbol);
        if (it != vdw_radii_.end()) {
            vdw_radius = it->second;
        } else {
            std::cout << "Warning: Unknown element " << symbol << ", using default radius 1.7"
                      << std::endl;
        }
        std::string name = symbol + std::to_string(i);
        out << " " << pos.x << " " << pos.y << " " << pos.z << " " << vdw_radius << " "
            << "1"
            << " " << name << std::endl;
    }

    out.close();
}

void Surface::savePointCloud(
    const std::string& filename, bool include_normals, bool include_colors,
    const std::vector<std::string>& properties,
    const std::unordered_map<std::string, std::any>* property_cache) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    // PLY header
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "comment Generated by SPALT molecular surface export\n";
    file << "element vertex " << vertices.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";

    if (include_normals) {
        file << "property float nx\n";
        file << "property float ny\n";
        file << "property float nz\n";
    }

    if (include_colors && !properties.empty()) {
        // Add a color column for each property
        for (const auto& property : properties) {
            file << "property float " << property << "_color\n";
        }
    }

    file << "end_header\n";

    for (const auto& [vertex_index, vertex_data] : vertices) {
        const auto& [atom_name, coords_normals] = vertex_data;

        file << coords_normals[0] << " " << coords_normals[1] << " " << coords_normals[2];

        if (include_normals) {
            file << " " << coords_normals[3] << " " << coords_normals[4] << " "
                 << coords_normals[5];
        }

        if (include_colors && !properties.empty()) {
            // Write color value for each property
            for (const auto& property : properties) {
                double color_value =
                    getVertexColorValue(property, vertex_index, atom_name, property_cache);
                file << " " << color_value;
            }
        }

        file << "\n";
    }

    file.close();
}

void Surface::saveMesh(const std::string& filename, bool include_normals, bool include_colors,
                       const std::vector<std::string>& properties,
                       const std::unordered_map<std::string, std::any>* property_cache) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    // PLY header
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "comment Generated by SPALT molecular surface export\n";
    file << "element vertex " << vertices.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";

    if (include_normals) {
        file << "property float nx\n";
        file << "property float ny\n";
        file << "property float nz\n";
    }

    if (include_colors && !properties.empty()) {
        // Add a color column for each property
        for (const auto& property : properties) {
            file << "property float " << property << "_color\n";
        }
    }

    file << "element face " << faces.size() << "\n";
    file << "property list uchar int vertex_indices\n";

    file << "end_header\n";

    for (const auto& [vertex_index, vertex_data] : vertices) {
        const auto& [atom_name, coords_normals] = vertex_data;

        file << coords_normals[0] << " " << coords_normals[1] << " " << coords_normals[2];

        if (include_normals) {
            file << " " << coords_normals[3] << " " << coords_normals[4] << " "
                 << coords_normals[5];
        }

        if (include_colors && !properties.empty()) {
            // Write color value for each property
            for (const auto& property : properties) {
                double color_value =
                    getVertexColorValue(property, vertex_index, atom_name, property_cache);
                file << " " << color_value;
            }
        }

        file << "\n";
    }

    for (const auto& face : faces) {
        file << face.size();
        for (int vertex_index : face) {
            file << " " << vertex_index;
        }
        file << "\n";
    }

    file.close();
}

void Surface::transformSurface(const Eigen::Isometry3d& transformation) {
    // Transform all vertices
    for (auto& vertex_pair : vertices) {
        auto& coords = vertex_pair.second.second;
        if (coords.size() >= 6) {  // x, y, z, nx, ny, nz
            // Transform position
            Eigen::Vector3d position(coords[0], coords[1], coords[2]);
            Eigen::Vector3d transformed_position = transformation * position;
            coords[0] = transformed_position.x();
            coords[1] = transformed_position.y();
            coords[2] = transformed_position.z();

            // Transform normal (rotation only, no translation)
            Eigen::Vector3d normal(coords[3], coords[4], coords[5]);
            Eigen::Vector3d transformed_normal = transformation.rotation() * normal;
            coords[3] = transformed_normal.x();
            coords[4] = transformed_normal.y();
            coords[5] = transformed_normal.z();
        }
    }
}

void Surface::generateFibonacciSurface(const spalt::FibonacciSurfaceParams& fibonacci_params) {
    // Use the EXACT same pattern as working Fibonacci tests
    auto mol_ptr = molecule.get_mol();
    if (!mol_ptr) {
        throw std::runtime_error("Molecule is null");
    }
    
    if (mol_ptr->getNumConformers() == 0) {
        throw std::runtime_error("Molecule has no conformers for surface generation");
    }
    
    const auto& conformer = mol_ptr->getConformer(conformer_id);
    int num_atoms = conformer.getNumAtoms();
    
    if (num_atoms == 0) {
        throw std::runtime_error("Conformer has no atoms");
    }
    
    // Extract atom centers and radii - EXACTLY like working test
    Eigen::MatrixXd centers(num_atoms, 3);
    Eigen::VectorXd radii(num_atoms);
    
    for (int i = 0; i < num_atoms; ++i) {
        const auto& pos = conformer.getAtomPos(i);
        centers(i, 0) = pos.x;
        centers(i, 1) = pos.y;
        centers(i, 2) = pos.z;
        
        // Get atom symbol and look up radius
        const auto& atom = mol_ptr->getAtomWithIdx(i);
        std::string symbol = atom->getSymbol();
        auto it = vdw_radii_.find(symbol);
        radii(i) = (it != vdw_radii_.end()) ? it->second : DEFAULT_VDW_RADIUS;
    }
    
    // Create Fibonacci surface generator
    spalt::FibonacciSurfaceGenerator generator(fibonacci_params);
    
    // Generate surface
    std::shared_ptr<spalt::SimpleMesh> mesh;
    std::shared_ptr<spalt::SimplePointCloud> point_cloud;
    auto result = generator.generateSurface(centers, radii);
    mesh = result.first;
    point_cloud = result.second;
    
    if (!mesh || !point_cloud) {
        std::cerr << "Error: Failed to generate Fibonacci surface - mesh or point_cloud is null" << std::endl;
        throw std::runtime_error("Failed to generate Fibonacci surface");
    }    
    
    // Convert simple mesh to internal format
    vertices.clear();
    faces.clear();
    
    // Convert vertices
    const auto& mesh_vertices = mesh->vertices;
    const auto& mesh_faces = mesh->triangles;
    
    for (size_t i = 0; i < mesh_vertices.size(); ++i) {
        const auto& vertex = mesh_vertices[i];
        
        // Compute normal from mesh faces if available
        Eigen::Vector3d normal(0.0, 0.0, 0.0);
        if (!mesh_faces.empty()) {
            normal = computeVertexNormal(i, mesh_vertices, mesh_faces);
        }
        
        std::vector<double> coords_normals = {
            vertex.x(), vertex.y(), vertex.z(), normal.x(), normal.y(), normal.z()
        };
        
        // Find closest atom for naming
        int closest_atom = 0;
        double min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < num_atoms; ++j) {
            double dist = (centers.row(j).transpose() - vertex).norm();
            if (dist < min_dist) {
                min_dist = dist;
                closest_atom = j;
            }
        }
        
        const auto& atom = molecule.get_mol()->getAtomWithIdx(closest_atom);
        std::string atom_name = atom->getSymbol() + std::to_string(closest_atom);
        
        vertices[i] = std::make_pair(atom_name, coords_normals);
    }
    
    // Convert faces (triangles)
    for (const auto& face : mesh_faces) {
        std::vector<int> face_indices = {face.x(), face.y(), face.z()};
        faces.push_back(face_indices);
    }
    
    // Apply subsampling if needed
    if (sample_method != "full") {
        applySubsampling();
    }
}

Eigen::Vector3d Surface::computeVertexNormal(size_t vertex_index, 
                                            const std::vector<Eigen::Vector3d>& vertices, 
                                            const std::vector<Eigen::Vector3i>& faces) const {
    Eigen::Vector3d normal_sum(0.0, 0.0, 0.0);
    int face_count = 0;
    
    // Safety check: ensure vertex_index is within bounds
    if (vertex_index >= vertices.size()) {
        std::cerr << "Warning: Invalid vertex_index " << vertex_index 
                  << " (total vertices: " << vertices.size() << ")" << std::endl;
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    
    // Find all faces that contain this vertex
    for (const auto& face : faces) {
        // Safety check: ensure face indices are within bounds
        if (face.x() < 0 || static_cast<size_t>(face.x()) >= vertices.size() ||
            face.y() < 0 || static_cast<size_t>(face.y()) >= vertices.size() ||
            face.z() < 0 || static_cast<size_t>(face.z()) >= vertices.size()) {
            continue;  // Skip invalid face
        }
        
        if (static_cast<size_t>(face.x()) == vertex_index || static_cast<size_t>(face.y()) == vertex_index || static_cast<size_t>(face.z()) == vertex_index) {
            // Compute face normal
            Eigen::Vector3d v0 = vertices[face.x()];
            Eigen::Vector3d v1 = vertices[face.y()];
            Eigen::Vector3d v2 = vertices[face.z()];
            
            Eigen::Vector3d edge1 = v1 - v0;
            Eigen::Vector3d edge2 = v2 - v0;
            Eigen::Vector3d face_normal = edge1.cross(edge2);
            
            // Normalize face normal
            double norm = face_normal.norm();
            if (norm > 1e-8) {
                face_normal /= norm;
                normal_sum += face_normal;
                face_count++;
            }
        }
    }
    
    // Average the face normals
    if (face_count > 0) {
        normal_sum /= face_count;
        // Normalize the result
        double norm = normal_sum.norm();
        if (norm > 1e-8) {
            normal_sum /= norm;
        }
    }
    
    return normal_sum;
}
