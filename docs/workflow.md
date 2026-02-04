# Workflow Guide

Complete workflow for generating and analyzing twisted layered material structures.

---

## Overview

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  1. Structure   │────▶│  2. LAMMPS      │────▶│  3. Analysis    │
│    Generation   │     │    Simulation   │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

---

## Step 1: Structure Generation

### 1.1 Prepare Input Files

Create the following files in your working directory:

```
MySystem/
├── twisterase.inp    # Main configuration
├── layer1.inp        # Bottom layer
├── layer2.inp        # Second layer
├── ...
└── potential files   # .sw, .KC, .tersoff, etc.
```

### 1.2 Run TwisterASE

```bash
cd Test/MySystem
python ../../Src/twisterase.py
```

### 1.3 Output Files

| File | Description |
|------|-------------|
| `superlattice.cif` | Full twisted structure (for visualization) |
| `superlattice_lammps.cif` | LAMMPS-compatible CIF |
| `structure.lammps` | LAMMPS data file |
| `lammps.in` | LAMMPS input script |

---

## Step 2: LAMMPS Simulation

### 2.1 Run Relaxation

```bash
lmp_serial -in lammps.in
```

Or for parallel execution:
```bash
mpirun -np 4 lmp_mpi -in lammps.in
```

### 2.2 Output Files

| File | Description |
|------|-------------|
| `dump.minimization` | Trajectory during minimization |
| `dump.Final` | Final relaxed structure |
| `log.lammps` | LAMMPS log file |

### 2.3 Required LAMMPS Packages

Ensure your LAMMPS build includes:
- `MANYBODY` (for sw/mod, tersoff)
- `INTERLAYER` (for kolmogorov/crespi/z)
- `MOLECULE` (for full atom style)

---

## Step 3: Post-Processing

### 3.1 Create cutpos.inp

```python
n_layers = 2
lammps_dump = "dump.Final"
```

### 3.2 Extract Layers

```bash
python ../../Src/cutpos.py
```

**Output:**
- `relaxed_structure.cif` - Full relaxed structure
- `layer_1.cif`, `layer_2.cif`, ... - Individual layers
- `layer_1_coords.dat`, `layer_2_coords.dat`, ... - Coordinate files

---

## Step 4: Analysis

### 4.1 Automated Analysis

```bash
python ../../Src/run_analysis.py
```

This runs:
1. Directory setup
2. Input file generation
3. Interlayer spacing analysis
4. Strain analysis

### 4.2 Output Structure

```
MySystem/
├── InterlayerSpacingMap/
│   └── Layer_1-2/
│       ├── scatter.png
│       └── interpolated.png
└── StrainMap/
    ├── Layer_1/
    │   ├── strain_Layer_1.png
    │   └── hist_strain_Layer_1.png
    └── Layer_2/
        ├── strain_Layer_2.png
        └── hist_strain_Layer_2.png
```

---

## Step 5: DFT Conversion (Optional)

### SIESTA Format

```bash
python ../../Src/cif2siesta.py layer_1.cif layer_1.fdf
```

### Quantum ESPRESSO Format

```bash
python ../../Src/cif2qe.py layer_1.cif layer_1.in
```

---

## Quick Reference

### Complete Workflow Commands

```bash
# 1. Generate structure
python ../../Src/twisterase.py

# 2. Run LAMMPS
lmp_serial -in lammps.in

# 3. Extract layers
python ../../Src/cutpos.py

# 4. Run analysis
python ../../Src/run_analysis.py
```

### Troubleshooting

| Issue | Solution |
|-------|----------|
| LAMMPS potential error | Check potential files are in working directory |
| Missing atom types | Verify tags are unique across all layers |
| Wrong interlayer spacing | Check translate_z values in layer files |
| Structure not periodic | Verify superlattice vectors match twist angle |
