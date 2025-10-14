# Development Guide for spalt

This document provides comprehensive information for developers working on the spalt molecular surface generation library.

## 🏗️ Architecture Overview

### Core Design Principles

The library follows a **molecule-centric design** with complete encapsulation:

- **Molecule**: Primary entity for molecular structure representation, surface generation, and property computation
- **Surface**: Molecular surface representation with multi-property coloring capabilities
- **ISurfaceProperty**: Extensible property computation framework with generic access
- **PropertyParamsRegistry**: Centralized parameter management for consistent color mapping
- **Dynamic Property Discovery**: Automatic detection and integration of all available properties

### API Design Evolution

The library has evolved from a surface-centric design to a **molecule-centric design** with **complete encapsulation**:

**Original API (Surface-centric - unsafe):**
```cpp
auto molecule = std::make_shared<Molecule>("CCO");
Surface surface(molecule, 1000);  // Surface owns/uses molecule
surface->registerProperty(...);   // Direct surface manipulation
```

**Current API (Molecule-centric - safe):**
```cpp
Molecule molecule("CCO");                                    // Molecule is primary entity
auto property = molecule.compute("esp");                    // Direct property computation
auto esp_values = molecule.getProperty<std::vector<double>>("esp");  // Type-safe access
```

**Key Safety Benefits:**
- **No Raw Pointers**: Surfaces and properties are managed internally by the molecule
- **Automatic Lifecycle**: Resources are automatically cleaned up when molecule is destroyed
- **Complete Encapsulation**: All operations go through the molecule interface
- **Intelligent Caching**: Property results are cached and invalidated automatically
- **Exception Safety**: Proper error handling and resource cleanup throughout

## 🎨 Multi-Property Colored Point Cloud Alignment

The library features a sophisticated colored point cloud alignment system that uses **ALL computed surface properties simultaneously** for enhanced alignment accuracy.

### How It Works

1. **Dynamic Property Discovery**: The system automatically detects all computed properties from the property cache
2. **Multi-Property Color Synthesis**: Each surface point gets a combined color from ALL available properties:
   ```cpp
   // Each property contributes to the final color
   for (const auto& property : available_properties) {
       double normalized_value = normalizePropertyValue(property, values[i]);
       
       // Property-specific color signature based on name hash
       std::hash<std::string> hasher;
       size_t hash = hasher(property);
       
       // Add property contribution to combined color
       combined_color[0] += r_factor * (1.0 + normalized_value) / 2.0;
       combined_color[1] += g_factor * (1.0 + normalized_value) / 2.0;
       combined_color[2] += b_factor * (1.0 + normalized_value) / 2.0;
   }
   ```
3. **Colored ICP Registration**: Open3D's `RegistrationColoredICP` uses the multi-property colors for superior alignment
4. **Automatic Fallback**: Falls back to standard ICP if colored ICP fails

### Benefits

- **Enhanced Accuracy**: Multi-dimensional property information improves alignment quality
- **Automatic Integration**: No manual configuration - all properties are used automatically
- **Consistent Color Mapping**: All properties normalized to [-1, 1] range for balanced contribution
- **Robust Fallback**: Graceful degradation to standard ICP if needed

### Example Output
```
Available properties for alignment: esp hbond hydrophobicity 
Attempting colored ICP registration with 3 properties...
✓ Used colored ICP registration with property-based colors (fitness: 1.000)
```

## 🔧 Implementing New Surface Properties

The library provides a flexible framework for implementing new surface properties. Here's a complete guide:

### 1. Create Property Header File

Create a new header file in `src/props/` (e.g., `my_property.h`):

```cpp
#pragma once
#include "surface_property.h"
#include <vector>
#include <string>

/**
 * @class MyProperty
 * @brief Description of your property
 */
class MyProperty : public ISurfaceProperty {
public:
    MyProperty() = default;
    ~MyProperty() = default;
    
    /**
     * @brief Get the property key identifier
     * @return String identifier for your property
     */
    std::string key() const override { return "my_property"; }
    
    /**
     * @brief Compute the property values
     * @param surface Surface object containing the molecular structure
     * @param cache Property cache for storing computed results
     */
    void compute(const Surface& surface, std::unordered_map<std::string, std::any>& cache) override;
    
private:
    // Add any private helper methods here
};
```

### 2. Create Property Implementation File

Create the implementation file `src/props/my_property.cpp`:

```cpp
#include "my_property.h"
#include "surface.h"
#include "property_params.h"
#include <iostream>

void MyProperty::compute(const Surface& surface, std::unordered_map<std::string, std::any>& cache) {
    try {
        // Get the molecule
        auto mol = surface.molecule.get_mol();
        if (!mol) {
            throw std::runtime_error("Failed to get RDKit molecule");
        }
        
        // Get surface point coordinates
        std::vector<std::vector<double>> surface_points;
        for (const auto& [vertex_index, vertex_data] : surface.vertices) {
            const auto& coords_normals = vertex_data.second;
            surface_points.push_back({coords_normals[0], coords_normals[1], coords_normals[2]});
        }
        
        // Compute your property values
        std::vector<double> property_values;
        for (size_t i = 0; i < surface_points.size(); ++i) {
            // Your property computation logic here
            double value = computeMyProperty(mol, surface_points[i]);
            property_values.push_back(value);
        }
        
        // Store results in cache
        cache["my_property_values"] = property_values;
        
        // Compute statistics for display
        if (!property_values.empty()) {
            double min_val = *std::min_element(property_values.begin(), property_values.end());
            double max_val = *std::max_element(property_values.begin(), property_values.end());
            double avg_val = std::accumulate(property_values.begin(), property_values.end(), 0.0) / property_values.size();
            
            cache["my_property_avg"] = avg_val;
            
            std::cout << "My property computation successful: " << property_values.size() << " values" << std::endl;
            std::cout << "  Actual min value: " << min_val << std::endl;
            std::cout << "  Actual max value: " << max_val << std::endl;
            std::cout << "  Average value: " << avg_val << std::endl;
        }
        
        // Store standardized min/max values for consistent color mapping across molecules
        PropertyParams my_property_params = PropertyParamsRegistry::getParams("my_property");
        cache["my_property_min"] = my_property_params.min_value;
        cache["my_property_max"] = my_property_params.max_value;
        
    } catch (const std::exception& e) {
        std::cerr << "Error computing my property: " << e.what() << std::endl;
        throw;
    }
}

// Private helper method implementation
double MyProperty::computeMyProperty(std::shared_ptr<RDKit::ROMol> mol, const std::vector<double>& point) {
    // Implement your specific property calculation here
    // This is where you'd put the core logic for your property
    
    // Example: simple distance-based property
    double total_distance = 0.0;
    int atom_count = 0;
    
    if (mol->getNumConformers() > 0) {
        const RDKit::Conformer& conf = mol->getConformer(0);
        for (int i = 0; i < conf.getNumAtoms(); ++i) {
            const RDGeom::Point3D& atom_pos = conf.getAtomPos(i);
            double dx = point[0] - atom_pos.x;
            double dy = point[1] - atom_pos.y;
            double dz = point[2] - atom_pos.z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            total_distance += distance;
            atom_count++;
        }
    }
    
    return atom_count > 0 ? total_distance / atom_count : 0.0;
}
```

### 3. Register Property Parameters

Add your property to the parameter registry in `src/props/property_params.cpp`:

```cpp
// Add to the registry map
const std::unordered_map<std::string, PropertyParams> PropertyParamsRegistry::registry_ = {
    {"esp", PropertyParams("esp", -10.0, 10.0)},           // Electrostatic potential in eV
    {"hydrophobicity", PropertyParams("hydrophobicity", -4.5, 4.5)}, // Hydrophobicity scale
    {"hbond", PropertyParams("hbond", -1.0, 1.0)},         // Hydrogen bond potential
    {"my_property", PropertyParams("my_property", -5.0, 5.0)}  // Add your property here
};
```

**Important**: These min/max values are used for consistent color mapping across different molecules. Choose values that represent the typical range for your property.

### 4. Register Property in Molecule Class

Add your property to the molecule's property registry in `src/molecule.cpp`:

```cpp
void Molecule::registerDefaultProperties() {
    // Existing properties...
    registerProperty("esp", []() { 
        return std::make_unique<SurfaceESPProperty>(ChargeMethod::RDKIT); 
    });
    registerProperty("hydrophobicity", []() { 
        return std::make_unique<HydrophobicityProperty>(); 
    });
    registerProperty("hbond", []() { 
        return std::make_unique<HBondProperty>(); 
    });
    
    // Add your property
    registerProperty("my_property", []() { 
        return std::make_unique<MyProperty>(); 
    });
}
```

### 5. Update Build System

Add your source files to `CMakeLists.txt`:

```cmake
target_sources(spalt PRIVATE
  src/props/esp.cpp
  src/props/hydrophobicity.cpp
  src/props/hbond.cpp
  src/props/rdkit_utils.cpp
  src/props/property_params.cpp
  src/props/my_property.cpp  # Add your property here
  src/props/pharmacophore.cpp
)
```

And to `tests/CMakeLists.txt`:

```cmake
target_sources(spalt_tests PRIVATE
  # ... existing sources ...
  ${CMAKE_SOURCE_DIR}/src/props/my_property.cpp  # Add your property here
)
```

### 6. Add CLI Support (Optional)

To make your property available via the command line, update `src/cli.cpp`:

```cpp
void printUsage() {
    std::cout << "  --properties P    Surface properties to compute (default: esp,hydrophobicity,hbond):" << std::endl;
    std::cout << "                      - esp: Electrostatic potential" << std::endl;
    std::cout << "                      - hydrophobicity: Hydrophobicity values" << std::endl;
    std::cout << "                      - hbond: Hydrogen bond potential" << std::endl;
    std::cout << "                      - my_property: My custom property" << std::endl;  // Add your property
    std::cout << "                      - all: Compute all available properties" << std::endl;
}
```

### 7. Create Tests

Create test file `tests/test_my_property.cpp`:

```cpp
#include <gtest/gtest.h>
#include "molecule.h"
#include "my_property.h"

class MyPropertyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test molecule
        molecule = std::make_unique<Molecule>("CCO");  // Ethanol
        surface_id = molecule->createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms");
    }
    
    std::unique_ptr<Molecule> molecule;
    int surface_id;
};

TEST_F(MyPropertyTest, ComputeMyProperty) {
    // Test property computation
    molecule->compute("my_property", surface_id);
    
    // Verify property was computed
    auto values = molecule->getProperty<std::vector<double>>("my_property", "my_property_values", surface_id);
    EXPECT_FALSE(values.empty());
    
    // Test property values are reasonable
    for (double value : values) {
        EXPECT_TRUE(std::isfinite(value));
        EXPECT_GE(value, 0.0);  // Assuming positive values
    }
}

TEST_F(MyPropertyTest, PropertyCaching) {
    // Test that property values are cached
    molecule->compute("my_property", surface_id);
    auto values1 = molecule->getProperty<std::vector<double>>("my_property", "my_property_values", surface_id);
    
    molecule->compute("my_property", surface_id);
    auto values2 = molecule->getProperty<std::vector<double>>("my_property", "my_property_values", surface_id);
    
    EXPECT_EQ(values1.size(), values2.size());
    for (size_t i = 0; i < values1.size(); ++i) {
        EXPECT_DOUBLE_EQ(values1[i], values2[i]);
    }
}
```

### 8. Usage Examples

Once implemented, your property can be used like any other:

```cpp
// C++ API
Molecule molecule("CCO");
auto property = molecule.compute("my_property");
auto values = molecule.getProperty<std::vector<double>>("my_property", "my_property_values");

// CLI
./spalt reference.sdf input.sdf output/ --properties my_property
```

### Key Design Principles

1. **Inherit from ISurfaceProperty**: All properties must implement the base interface
2. **Use Property Cache**: Store computed values in the provided cache with consistent naming
3. **Standardize Min/Max**: Use PropertyParamsRegistry for consistent color mapping
4. **Error Handling**: Wrap computations in try-catch blocks with meaningful error messages
5. **Statistics**: Compute and display min/max/average values for debugging
6. **Testing**: Create comprehensive tests for your property
7. **Documentation**: Add proper Doxygen comments

### Property Cache Keys

Use these consistent naming patterns for cache keys:
- `{property_name}_values`: Main property values vector
- `{property_name}_min`: Minimum value for color mapping
- `{property_name}_max`: Maximum value for color mapping  
- `{property_name}_avg`: Average value for statistics

### Color Mapping

The library automatically handles color mapping using the min/max values you provide. Values are normalized to the [-1, 1] range for consistent visualization across all properties.

## 📊 PLY Output Format

The library exports molecular surfaces in PLY format with support for multiple color columns, one for each computed property.

### PLY Header Example
```
ply
format ascii 1.0
comment Created by spalt molecular surface generator
element vertex 1000
property float x
property float y
property float z
property float nx
property float ny
property float nz
property float esp_color
property float hydrophobicity_color
property float hbond_color
element face 2000
property list uchar uint vertex_indices
end_header
```

### Color Column Format
- Each property gets its own `float {property_name}_color` column
- Values range from -1.0 to 1.0
- -1.0 typically represents minimum property values (red in visualization)
- 1.0 typically represents maximum property values (blue in visualization)
- 0.0 represents neutral/average values

### Visualization in PyMOL
```python
# Load the surface
cmd.load("surface.ply")

# Color by specific property
cmd.spectrum("esp_color", "red_white_blue", "surface")

# Show multiple properties simultaneously
cmd.spectrum("hydrophobicity_color", "red_white_blue", "surface")
```

## 🧪 Testing

### Test Suite Overview
The project includes **69 comprehensive tests** covering all major functionality:

#### **Test Categories**
- **Surface Generation Tests** (1 test): Surface generation with different sampling methods
- **Utility Tests** (2 tests): File I/O and string processing utilities
- **Molecule Tests** (17 tests): Molecule construction, conformer generation, and basic operations
- **ESP Property Tests** (10 tests): Electrostatic potential computation and caching
- **Hydrogen Bond Tests** (4 tests): HBond property computation and feature detection
- **Hydrophobicity Tests** (2 tests): Hydrophobicity calculation and LogP conversion
- **Conformer Generation Tests** (10 tests): Advanced conformer generation with clustering
- **Property System Tests** (13 tests): Dynamic property system and parameter management
- **Alignment Tests** (10 tests): Multi-property colored alignment and transformation

#### **Test Coverage Areas**
✅ **Core Functionality**
- Surface generation with MSMS integration
- Molecule construction from SMILES and coordinate files
- Conformer generation (simple and advanced)
- Property computation (ESP, hydrophobicity, hbond)

✅ **Advanced Features**
- Advanced conformer generation with K-means clustering
- Multi-property colored alignment using Open3D
- Dynamic property registration and discovery
- Conformer-surface correspondence verification

✅ **Edge Cases & Error Handling**
- Invalid input handling
- Boundary condition testing
- Resource management and cleanup
- Exception handling and fallback mechanisms

#### **Running Tests**
```bash
# Build and run all tests
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
ctest --output-on-failure

# Run specific test categories
ctest -R "ConformerGenerationTest"
ctest -R "PropertySystemTest"
ctest -R "AlignmentTest"

# Run with verbose output
ctest --output-on-failure --verbose
```

#### **CI/CD Pipeline**
- **Automated Testing**: Runs on multiple compilers (GCC 9-11, Clang 12-13)
- **Build Matrix**: Tests Debug and Release builds
- **Code Quality**: Static analysis with clang-tidy and cppcheck
- **Coverage Reporting**: Automated coverage collection and reporting
- **Documentation**: Automated documentation generation
- **Security Scanning**: Vulnerability scanning with Trivy

### Test Results
```
100% tests passed, 0 tests failed out of 69
Total Test time (real) = 2.82 sec
```

## 🏗️ Project Structure
```
spalt/
├── src/                           # Source code
│   ├── molecule.h/cpp            # Main molecular wrapper class
│   ├── surface.h/cpp             # Surface generation and manipulation
│   ├── cli.h/cpp                 # Command-line interface
│   ├── utils.h/cpp               # Utility functions
│   ├── fps.h/cpp                 # Farthest Point Sampling algorithm
│   ├── conformer_generator.h/cpp # Advanced conformer generation
│   ├── surface_fibonacci.h/cpp   # Fibonacci sphere surface generation
│   └── props/                    # Surface properties
│       ├── surface_property.h    # Base property interface
│       ├── property_params.h/cpp # Property parameter registry
│       ├── esp.h/cpp             # Electrostatic Potential
│       ├── hydrophobicity.h/cpp  # Hydrophobicity (BRICS-based)
│       ├── hbond.h/cpp           # Hydrogen Bond Potential
│       ├── rdkit_utils.h/cpp     # RDKit feature utilities
│       └── pharmacophore.h/cpp   # Pharmacophore features (in development)
├── tests/                        # GoogleTest-based tests
│   ├── test_molecule.cpp         # Molecule functionality tests
│   ├── test_surface.cpp          # Surface generation tests
│   ├── test_esp.cpp              # ESP computation tests
│   ├── test_hydrophobicity.cpp   # Hydrophobicity tests
│   ├── test_hbond_features.cpp   # Hydrogen bond feature tests
│   ├── test_conformer_generation.cpp # Conformer generation tests
│   ├── test_property_system.cpp  # Property system tests
│   ├── test_alignment.cpp        # Alignment tests
│   ├── test_fibonacci_surface.cpp # Fibonacci surface tests
│   └── test_utils.cpp            # Utility function tests
├── examples/                     # Example input data
├── context/                      # Reference implementations and documentation
├── scripts/                      # Utility scripts
├── CMakeLists.txt                # Main build configuration
├── README.md                     # User-focused documentation
├── DEVELOPMENT.md                # Developer documentation (this file)
├── CONTRIBUTING.md               # Contributing guidelines
└── LICENSE                       # MIT License
```

## 🎯 TODO

### High Priority
- [ ] **Pharmacophore Implementation**: Complete the pharmacophore property computation
- [ ] **Performance Optimization**: Optimize surface generation and property computation for large molecules
- [ ] **Documentation Generation**: Set up automated Doxygen documentation generation
- [ ] **Package Management**: Add vcpkg/Conan package manager support

### Medium Priority
- [ ] **Export Formats**: Add more surface export formats (OBJ, STL, etc.) and conformer formats
- [ ] **Alignment Scoring**: Implement alignment scoring and quality metrics for surface alignment
- [ ] **Documentation**: Generate and maintain comprehensive API documentation
- [ ] **Batch Processing**: Add CLI options for processing multiple molecules in batch

### Low Priority
- [ ] **Additional Properties**: Implement more surface properties (shape, volume, etc.)
- [ ] **Parallel Processing**: Add multi-threading support for batch operations
- [ ] **Memory Optimization**: Optimize memory usage for large-scale molecular datasets
- [ ] **GUI Interface**: Consider adding a simple GUI for interactive molecular visualization

### Completed Features ✅
- [x] **Molecule-Centric API**: Complete refactor to molecule-centric design
- [x] **Dynamic Property System**: Fully extensible system with automatic property discovery and registration
- [x] **Multi-Property Colored Alignment**: Simultaneous use of ALL computed properties for enhanced point cloud alignment
- [x] **Property Parameter Registry**: Centralized parameter management for consistent color mapping across molecules
- [x] **ESP Computation**: XTB and RDKit-based electrostatic potential calculation
- [x] **Hydrophobicity Calculation**: BRICS fragmentation-based LogP calculation with literature fallback
- [x] **Hydrogen Bond Potential**: Complete donor/acceptor detection with geometric calculations
- [x] **Advanced Conformer Generation**: K-means clustering, MMFF optimization, and diverse sampling strategies
- [x] **Multi-Conformer Alignment**: Process and align multiple conformers with individual surface representations
- [x] **Surface Subsampling**: Random, FPS, and full surface sampling methods
- [x] **Advanced CLI Interface**: Command-line tool with flexible property selection, conformer control, and output options
- [x] **Multi-Color PLY Export**: Separate color columns for each property in PLY files (-1 to 1 range)
- [x] **SDF Export with Alignment Scores**: Embed fitness and RMSE scores directly into SDF files
- [x] **Text File SMILES Support**: Process multiple SMILES from text files with batch conformer generation
- [x] **Fibonacci Surface Generation**: Alternative surface generation method using Open3D's ball pivoting
- [x] **File I/O**: PLY surface export and SDF conformer export
- [x] **Comprehensive Testing**: 69 tests covering all major functionality with 100% pass rate
- [x] **No Hardcoded Values**: Fully dynamic system with no hardcoded property names or values
- [x] **CI/CD Pipeline**: Automated testing, coverage reporting, and code quality checks
- [x] **Professional Documentation**: Comprehensive README with badges and testing information

## 🔍 Debugging and Troubleshooting

### Common Issues

#### Build Issues
- **RDKit not found**: Ensure `-DRDKIT_ROOT` is set correctly for your installation
- **Open3D missing**: Install Open3D via conda or build from source
- **CMake version**: Requires CMake 3.16+ for proper C++17 support

#### Runtime Issues
- **Surface generation fails**: Check if MSMS is available or use Fibonacci method
- **Property computation errors**: Verify molecule has valid conformers
- **Alignment failures**: Ensure surfaces have sufficient vertices and properties are computed

#### Performance Issues
- **Slow surface generation**: Reduce vertex count or use simpler sampling methods
- **Memory issues**: Process molecules individually or reduce conformer count
- **Alignment timeouts**: Increase convergence criteria or reduce surface complexity

### Debug Output

Enable debug output by defining `DEBUG` macro:
```cpp
#define DEBUG
```

This will provide detailed information about:
- Property computation progress
- Alignment algorithm selection
- Surface generation statistics
- Caching behavior

### Profiling

Use standard profiling tools to identify bottlenecks:
```bash
# Profile with gprof
g++ -pg -O2 -g src/*.cpp
./spalt ...
gprof ./spalt gmon.out

# Profile with Valgrind
valgrind --tool=callgrind ./spalt ...
kcachegrind callgrind.out.*
```

## 📚 Additional Resources

### External Dependencies
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [Open3D Documentation](http://www.open3d.org/docs/)
- [Eigen Documentation](https://eigen.tuxfamily.org/dox/)
- [MSMS Documentation](http://mgltools.scripps.edu/packages/MSMS/)

### Related Projects
- [RDKit](https://github.com/rdkit/rdkit): Molecular informatics toolkit
- [Open3D](https://github.com/isl-org/Open3D): 3D data processing library
- [XTB](https://github.com/grimme-lab/xtb): Extended tight-binding program package

### Papers and References
- ESP calculation methods and validation
- BRICS fragmentation algorithm
- Hydrogen bond geometry calculations
- Point cloud registration algorithms

---

This development guide provides comprehensive information for extending and maintaining the spalt library. For user-focused documentation, see [README.md](README.md).


