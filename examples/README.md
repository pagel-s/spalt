# Example Input Files

This directory contains example molecular input files for testing and demonstrating the spalt library functionality.

## 📁 File Descriptions

### **SDF Files (Structure Data Format)**
- **`Mol1.sdf`** - Example molecule 1 in SDF format with 3D coordinates
- **`Mol2.sdf`** - Example molecule 2 in SDF format with 3D coordinates  
- **`Mol2_2d.sdf`** - Example molecule 2 in 2D SDF format (requires conformer generation)

### **SMILES Files**
- **`smiles.txt`** - Text file containing multiple SMILES strings (one per line)
  ```
  Cn1cc(-c2cnc(N)c3c(-c4ccc(Oc5ccccc5)cc4)csc23)cn1
  C1CCCCC1
  CCO
  c1ccccc1
  CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
  ```
- **`test_mol.smiles`** - Single SMILES string file for testing

### **XYZR Files (Coordinates + Radii)**
- **`test.xyzr`** - Example XYZR file with atomic coordinates and van der Waals radii
  - Format: `x y z radius element_name`
  - Used for direct surface generation without molecular structure parsing

### **Surface Files (Generated Output)**
- **`surface.vert`** - Example surface vertices file (13,489 vertices)
- **`surface.face`** - Example surface faces file (10,102 triangular faces)
  - These are example output files showing the surface generation format
  - Generated from molecular structures using MSMS or similar tools

## 🚀 Usage Examples

### **Using SDF Files**
```bash
# Generate surface from SDF file
./spalt examples/Mol1.sdf examples/Mol2.sdf output/ --properties esp,hydrophobicity

# Process 2D SDF (will generate conformers)
./spalt examples/Mol2_2d.sdf output/ --use-advanced --total-confs 20 --num-clusters 3
```

### **Using SMILES Files**
```bash
# Process single SMILES
./spalt examples/test_mol.smiles output/ --properties all

# Process multiple SMILES from text file
./spalt examples/Mol1.sdf examples/smiles.txt output/ --use-advanced --conformers 5
```

### **Using XYZR Files**
```bash
# Direct surface generation from XYZR
./spalt examples/test.xyzr output/ --properties esp
```

## 📊 File Formats

### **SDF Format**
- Standard molecular structure format
- Contains 3D coordinates, bonds, and properties
- Supported by most molecular visualization software

### **SMILES Format**
- Simplified molecular input line entry system
- Compact text representation of molecular structure
- Automatically converted to 3D coordinates by spalt

### **XYZR Format**
- Custom format: `x y z radius element_name`
- Direct atomic coordinates with van der Waals radii
- Bypasses molecular structure parsing for direct surface generation

## 🔬 Example Molecules

- **Mol1**: Complex heterocyclic compound with multiple rings
- **Mol2**: Aromatic compound with functional groups
- **SMILES examples**: Range from simple (benzene) to complex (drug-like molecules)

## 📝 Notes

- All example files are small and suitable for testing
- SDF files contain 3D coordinates ready for surface generation
- SMILES files will trigger conformer generation
- XYZR files bypass molecular parsing for direct surface generation
- Files are compatible with the spalt CLI interface
