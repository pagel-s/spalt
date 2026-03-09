# spalt

[![CI/CD Pipeline](https://github.com/username/spalt/actions/workflows/ci.yml/badge.svg)](https://github.com/username/spalt/actions/workflows/ci.yml)
[![Code Coverage](https://codecov.io/gh/username/spalt/branch/main/graph/badge.svg)](https://codecov.io/gh/username/spalt)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![RDKit](https://img.shields.io/badge/RDKit-Required-orange.svg)](https://www.rdkit.org/)
[![Open3D](https://img.shields.io/badge/Open3D-Required-red.svg)](http://www.open3d.org/)

**Molecular Surface Generation and Property Computation Library**

A C++ library for generating molecular surfaces from SMILES strings or coordinate files, with electrostatic potential calculations, surface property computation, and colored point cloud alignment.

<p align="center">
  <img src="assets/spalt.png" width="70%">
</p>

## Features

- **Molecular Surface Generation**: MSMS and Fibonacci sphere surface generation
- **Surface Properties**: Electrostatic potential (ESP), hydrophobicity, hydrogen bond potential
- **Advanced Conformer Generation**: K-means clustering and MMFF optimization
- **Multi-Conformer Alignment**: Process multiple conformers with individual surface representations
- **Configurable Charge Methods**: XTB (quantum chemistry) or RDKit (fast) charges for ESP
- **CLI Interface**: Command-line tool with flexible property selection and output options

## Requirements

- **CMake 3.16+**
- **C++17+ compiler** (GCC 9+, Clang 10+, MSVC 2019+)
- **RDKit (C++)**: Required for molecular structure handling
- **Eigen3**: Required for linear algebra operations
- **Open3D**: Required for point cloud alignment
- **MSMS**: Optional for surface generation (auto-detected)
- **XTB**: Optional for quantum chemistry calculations (auto-detected)

## Installation

### Prerequisites

**Important**: RDKit must be installed with C++ development headers. `pip install rdkit` does NOT include these headers.


### System Packages (Ubuntu/Debian)
```bash
# Install system dependencies (includes RDKit C++ headers)
sudo apt update
sudo apt install build-essential cmake ninja-build \
    librdkit-dev rdkit-data libeigen3-dev libopen3d-dev

# Clone and build
git clone https://github.com/username/spalt.git
cd spalt
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=17 \
  -DSPALT_BUILD_TESTS=ON \
  [-DRDKIT_ROOT=/path/to/rdkit/installation] # if cant be found/ manually installed

cmake --build build -j$(nproc)
cmake --build build -j$(nproc)
cmake --install build  # optional, installs library and executable
```

The CMake build system automatically detects RDKit installations from:
- Conda environments (`$CONDA_PREFIX`)
- Standard system locations (`/usr`, `/usr/local`, `/opt/rdkit`)
- Custom paths via `RDKIT_ROOT`

No manual configuration is needed for most installations!


## Example Usage

### Basic Alignment
```bash
# Basic alignment with default properties (all properties)
./spalt reference.sdf input.sdf output_dir/

# Align SMILES to reference
./spalt reference.sdf "CCO" output_dir/ --properties esp,hb,hy

# Use abbreviated property names
./spalt reference.sdf input.sdf output_dir/ --properties esp,hb,hy
```

### Advanced Conformer Generation
```bash
# Generate 10 diverse conformers from 100 samples and align best ones
./spalt reference.sdf "CCO" output_dir/ --conformers 10 --sample 100 --top-n 5

# Use XTB quantum chemistry charges for high-accuracy ESP
./spalt reference.sdf input.sdf output_dir/ --charge-method xtb --properties esp

# Use Fibonacci surface generation
./spalt reference.sdf input.sdf output_dir/ --mesh fibonacci --properties esp
```

### Output Options
```bash
# Export as point cloud with normals
./spalt reference.sdf input.sdf output_dir/ --output-type points --include-normals

# Save individual surface files
./spalt reference.sdf input.sdf output_dir/ --save-meshes

# Custom surface parameters
./spalt reference.sdf input.sdf output_dir/ --radius 1.5 --properties esp,hb
```

## CLI Options

### **Input/Output**
- `reference_file`: Reference molecule file (.sdf, .mol, .mol2, .pdb)
- `input`: Input molecule (file, SMILES string, or text file with SMILES). Text files can optionally contain a molecule name in the second column separated by a tab or space.
- `output_dir`: Output directory for results

### **Surface Properties**
- `--properties P`: Select properties (`esp`, `hb`, `hy`, `all`, `none`)
- `--charge-method M`: ESP charge method (`xtb`, `rdkit`) (default: rdkit)
- `--mesh {msms|fibonacci}`: Surface generation method (default: fibonacci)
- `--vertices N`: Number of surface vertices (default: 1000)
- `--radius R`: Probe radius for surface generation (default: 1.4)

### **Conformer Generation**
- `--conformers N`: Final number of conformers to align (default: 1)
- `--sample M`: Number of candidates to sample from (default: 1). If M > N, clustering is used to pick diverse conformers.
- `--use-advanced`: Use advanced conformer generation (shortcut for `--sample 50`)
- `--random-seed N`: Random seed for reproducible generation (default: 4)
- `--top-n N`: Save only the top N best aligned conformers (default: all)

### **Output Options**
- `--output-type {mesh|points}`: Export as triangular mesh or point cloud (default: mesh)
- `--include-normals`: Include vertex normals in output (default: false)
- `--save-meshes`: Save individual surface files for each conformer

## Testing

```bash
# Build and run all tests
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
ctest --output-on-failure
```


## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Add tests for new functionality
4. Ensure all tests pass (`ctest --output-on-failure`)
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

This does not include dependencies. For all dependencies for which a commercial license
is necessary, and alternative should be implemented.
 

## Acknowledgments

- **RDKit**: For molecular informatics capabilities
- **Open3D**: For 3D data processing and point cloud alignment
- **Eigen**: For linear algebra operations
- **MSMS**: For molecular surface generation
- **XTB**: For quantum chemistry calculations