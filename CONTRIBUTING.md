# Contributing to spalt

Thank you for your interest in contributing to spalt! This document provides guidelines and information for contributors.

## 🚀 Getting Started

### Prerequisites
- CMake 3.16+
- C++17 compatible compiler
- RDKit (C++)
- Eigen3
- Open3D

### Development Setup
```bash
# Clone the repository
git clone https://github.com/username/spalt.git
cd spalt

# Create build directory
mkdir build && cd build

# Configure with development settings
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic"

# Build
make -j$(nproc)

# Run tests
ctest --output-on-failure
```

## 📝 Code Style

We use Google C++ style with some modifications:
- Line length: 100 characters
- Indentation: 4 spaces
- Pointer alignment: Left (`int* ptr` not `int *ptr`)

### Formatting
```bash
# Format code with clang-format
find src/ tests/ -name "*.cpp" -o -name "*.h" | xargs clang-format -i
```

## 🧪 Testing

### Test Requirements
- All new features must include comprehensive tests
- Maintain 100% test pass rate
- Add tests for edge cases and error conditions

### Running Tests
```bash
# Run all tests
ctest --output-on-failure

# Run specific test categories
ctest -R "ConformerGenerationTest"
ctest -R "PropertySystemTest"

# Run with verbose output
ctest --output-on-failure --verbose
```

### Test Structure
- Tests should be in the `tests/` directory
- Use Google Test framework
- Follow naming convention: `ClassName.TestName`
- Use descriptive test names that explain what is being tested

## 🔧 Architecture Guidelines

For detailed architecture information, see [DEVELOPMENT.md](DEVELOPMENT.md).

### Key Principles
- **Molecule-Centric Design**: All operations should be performed through the `Molecule` class
- **Property System**: Use the extensible property system for new functionality
- **Error Handling**: Implement proper error handling and fallback mechanisms
- **Testing**: Maintain 100% test pass rate

## 📋 Pull Request Process

1. **Fork and Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Development**
   - Write code following the style guidelines
   - Add comprehensive tests
   - Update documentation as needed

3. **Testing**
   ```bash
   # Ensure all tests pass
   ctest --output-on-failure
   
   # Check code formatting
   find src/ tests/ -name "*.cpp" -o -name "*.h" | xargs clang-format --dry-run --Werror
   ```

4. **Commit**
   ```bash
   git add .
   git commit -m "Add feature: brief description"
   ```

5. **Submit PR**
   - Push to your fork
   - Create pull request with clear description
   - Link any related issues

### PR Requirements
- [ ] All tests pass
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Clear commit messages
- [ ] PR description explains changes

## 🐛 Bug Reports

When reporting bugs, please include:
- Clear description of the issue
- Steps to reproduce
- Expected vs actual behavior
- System information (OS, compiler, dependencies)
- Relevant code snippets or error messages

## 💡 Feature Requests

For new features:
- Check existing issues first
- Provide clear use case and motivation
- Consider implementation complexity
- Discuss in issues before starting work

## 📚 Documentation

### Code Documentation
- Use Doxygen-style comments for public APIs
- Document complex algorithms and design decisions
- Keep comments up-to-date with code changes

### README Updates
- Update feature lists when adding new functionality
- Keep installation instructions current
- Update examples and usage patterns

## 🔍 Code Review Process

### Reviewers
- All PRs require at least one review
- Focus on correctness, performance, and maintainability
- Check test coverage and edge cases

### Review Checklist
- [ ] Code follows project style
- [ ] Tests are comprehensive
- [ ] Documentation is updated
- [ ] No performance regressions
- [ ] Error handling is appropriate
- [ ] Memory management is correct

## 🏷️ Release Process

### Versioning
We use semantic versioning (MAJOR.MINOR.PATCH):
- MAJOR: Breaking API changes
- MINOR: New features, backward compatible
- PATCH: Bug fixes, backward compatible

### Release Checklist
- [ ] All tests pass
- [ ] Documentation updated
- [ ] Version numbers updated
- [ ] CHANGELOG.md updated
- [ ] Release notes prepared

## 🤝 Community Guidelines

### Code of Conduct
- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow
- Maintain professional communication

### Getting Help
- Check existing issues and documentation
- Ask questions in discussions
- Be specific about problems
- Provide context and examples

## 📞 Contact

- Create an issue for bugs or feature requests
- Use discussions for questions and ideas
- Follow the project for updates

Thank you for contributing to spalt! 🎉
