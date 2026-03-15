# SPALT

[![CI/CD Pipeline](https://github.com/pagel-s/spalt/actions/workflows/ci.yml/badge.svg)](https://github.com/pagel-s/spalt/actions/workflows/ci.yml)
[![Code Coverage](https://codecov.io/gh/pagel-s/spalt/branch/main/graph/badge.svg)](https://codecov.io/gh/pagel-s/spalt)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)

**Molecular Surface Generation and Property Computation Library**

A C++ library for generating molecular surfaces from SMILES strings or coordinate files, with electrostatic potential calculations, surface property computation, and colored point cloud alignment.

<p align="center">
  <img src="assets/spalt.png" width="70%" alt="SPALT Alignment Example">
</p>

## Features
- **Molecular Surface Generation**: High-performance Fibonacci sphere and MSMS surface generation.
- **Surface Properties**: Compute electrostatic potential (ESP), hydrophobicity, hydrogen bond potential, and 3D pharmacophore features.
- **Advanced Conformer Generation**: K-means clustering and MMFF optimization.
- **Multi-Conformer Alignment**: Process multiple conformers with individual surface representations concurrently via OpenMP.
- **Configurable Charges**: Support for XTB (quantum chemistry) or RDKit (fast) charges for ESP.
- **CLI Interface**: Powerful command-line tool with flexible property selection and mesh/point-cloud export.

## Requirements
- **CMake** 3.16+
- **C++17 Compiler** (GCC 9+, Clang 10+, MSVC 2019+)
- **RDKit (C++)**: Molecular structure handling
- **Eigen3**: Linear algebra operations
- **Open3D**: 3D data processing and point cloud alignment
- **OpenMP**: Highly recommended for parallel processing
- **Optional**: MSMS (surface generation), XTB (quantum chemistry charges)

> **Important**: RDKit must be installed with C++ development headers. `pip install rdkit` does NOT include these headers. Use `apt`, `conda`, or build from source.

## Building from Source

### 1. Install Dependencies (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install build-essential cmake ninja-build \
    librdkit-dev rdkit-data libeigen3-dev libopen3d-dev
```

### 2. Compile
```bash
git clone https://github.com/pagel-s/spalt.git
cd spalt

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# The binary is now available at build/spalt
./spalt --help
```
*Note: The CMake build system automatically detects RDKit installations from Conda environments, standard system locations, or via a custom `-DRDKIT_ROOT` flag.*

## Usage Examples

### Basic Alignment
```bash
# Basic alignment computing all surface properties
./spalt reference.sdf input.sdf output_dir/

# Align a SMILES string to the reference using abbreviated property names
./spalt reference.sdf "CCO" output_dir/ --properties esp,hb,hy
```

### Conformer Generation & Filtering
When providing 1D/2D input structures (e.g., SMILES strings or 2D SDF files), `spalt` dynamically generates 3D conformers before computing surfaces:
1. **Sampling**: Generate candidate conformers (`--sample`).
2. **Selection**: Cluster candidates and pick geometrically diverse structures (`--conformers`).
3. **Filtering**: Align and retain only the best-scoring structures (`--top-n`).

```bash
# Sample 100 conformers, select 10 diverse ones, align, and keep the top 5
./spalt reference.sdf "CCO" output_dir/ --sample 100 --conformers 10 --top-n 5

# Use a shortcut for robust sampling (--use-advanced is equivalent to --sample 50)
./spalt reference.sdf input.sdf output_dir/ --use-advanced --conformers 5

# Process multiple conformers using 4 CPU threads explicitly
./spalt reference.sdf "CCO" output_dir/ --use-advanced --conformers 10 --threads 4
```

### Output & Advanced Configs
```bash
# Export as point cloud with normal vectors
./spalt reference.sdf input.sdf output_dir/ --output-type points --include-normals

# Use XTB quantum chemistry charges for high-accuracy ESP
./spalt reference.sdf input.sdf output_dir/ --charge-method xtb --properties esp

# Output individual mesh files for each conformer processed
./spalt reference.sdf input.sdf output_dir/ --save-meshes
```

## CLI Reference

### Input/Output
- `reference`: Reference molecule (`.sdf`, `.mol`, `.mol2`, `.pdb`, or SMILES).
- `input`: Target molecule(s) (file, directory, SMILES string, or `.txt` with SMILES list).
- `output_dir`: Directory to save results.

### Surface Properties
- `--properties P`: Comma-separated list (`esp`, `hb`, `hy`, `pharma`, `all`, `none`).
- `--charge-method M`: ESP charge method (`xtb` or `rdkit`) [default: rdkit].
- `--mesh M`: Surface generation method (`msms` or `fibonacci`) [default: fibonacci].
- `--vertices N`: Number of surface vertices [default: 1000].
- `--radius R`: Probe radius in Ångstroms [default: 1.2].
- `--density D` / `--hdensity H`: MSMS surface density configurations [default: 3.0].
- `--type T`: Surface type (`tses` or `ases`) [default: tses].
- `--sample-method S`: Vertex sampling method (`full`, `fps`, `random`) [default: full].

### Conformer Pipeline
- `--sample M`: Number of candidates to sample and optimize [default: 1].
- `--conformers N`: Final number of diverse conformers to process [default: 1].
- `--top-n K`: Keep only the `K` best-aligned conformers [default: all].
- `--use-advanced`: Sets `--sample 50` for deeper conformational space searching.
- `--random-seed N`: Seed for reproducible generation [default: 4].

### Processing & Output
- `--threads N`: CPU threads for parallel multi-conformer processing [default: auto].
- `--addH` / `--removeH`: Control explicit hydrogen addition/removal.
- `--output-type T`: Output format (`mesh` or `points`) [default: mesh].
- `--include-normals`: Include vertex normals in outputs [default: false].
- `--save-meshes`: Save individual surface representations for each conformer.

## Testing
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DSPALT_BUILD_TESTS=ON
make -j$(nproc)
ctest --output-on-failure
```

## Contributing
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/new-feature`).
3. Add tests for new functionality.
4. Ensure all tests pass (`ctest --output-on-failure`).
5. Submit a pull request.

## License
Licensed under the [MIT License](LICENSE).
*Note: This does not include dependencies. For dependencies requiring a commercial license, an alternative should be implemented.*