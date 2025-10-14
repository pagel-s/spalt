# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Advanced conformer generation with K-means clustering and MMFF optimization
- Multi-conformer alignment with individual surface representations
- Text file SMILES support for batch processing
- SDF export with embedded alignment scores (fitness, RMSE)
- Comprehensive test suite with 69 tests covering all functionality
- Professional CI/CD pipeline with automated testing and coverage
- GitHub Actions workflows for CI, coverage, and code quality
- Professional README with badges and comprehensive documentation

### Changed
- Enhanced CLI interface with conformer control options
- Improved property parameter system with centralized management
- Updated alignment system to use all computed properties simultaneously
- Refactored conformer generation to be more configurable and robust

### Fixed
- Conformer-surface correspondence issues in multi-conformer workflows
- Property key naming inconsistencies in test suite
- Advanced conformer generation for text file inputs
- Edge case handling in conformer clustering

## [0.3.0] - 2024-01-XX

### Added
- Multi-property colored alignment using Open3D's Colored ICP
- Dynamic property system with automatic discovery and registration
- Property parameter registry for consistent color mapping
- Multi-color PLY export with separate columns for each property
- Advanced CLI interface with flexible property selection
- Hydrogen bond potential computation with geometric calculations
- BRICS fragmentation-based hydrophobicity calculation

### Changed
- Refactored to molecule-centric API design
- Improved surface caching with parameter-based invalidation
- Enhanced error handling and fallback mechanisms
- Updated property computation workflow

### Fixed
- Surface generation pipeline issues
- Property computation ordering problems
- Memory management and resource cleanup

## [0.2.0] - 2024-01-XX

### Added
- ESP computation with XTB and RDKit support
- Surface subsampling methods (random, FPS, full)
- PLY export for surfaces and meshes
- Basic conformer generation and alignment
- File I/O support for multiple formats

### Changed
- Improved surface generation reliability
- Enhanced coordinate file parsing
- Better error handling and validation

## [0.1.0] - 2024-01-XX

### Added
- Initial molecular surface generation
- MSMS integration with automatic detection
- Basic SMILES to 3D coordinate conversion
- Van der Waals radii database
- CMake build system
- Basic CLI interface

---

## Legend

- **Added** for new features
- **Changed** for changes in existing functionality
- **Deprecated** for soon-to-be removed features
- **Removed** for now removed features
- **Fixed** for any bug fixes
- **Security** for vulnerability fixes
