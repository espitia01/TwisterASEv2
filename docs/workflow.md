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
├── layer2.inp        # Second layer (add more as needed)
└── potential files   # .rebo, .KC, .sw, .tersoff, etc.
```

### 1.2 Run TwisterASE

```bash
cd Examples/MySystem
python ../../Src/twisterase.py
```

### 1.3 Output Files

| File | Description |
|------|-------------|
| `superlattice.cif` | Full twisted structure (for visualization) |
| `structure.lammps` | LAMMPS data file (`atom_style atomic`) |
| `lammps.in` | LAMMPS input script (ready to run) |
| `structure_ortho.lammps` | Orthorhombic LAMMPS data file (TMD ortho only) |

---

## Step 2: LAMMPS Simulation

### 2.1 Run Relaxation

```bash
lmp_serial -in lammps.in
```

Or for parallel execution:
```bash
mpirun -np 8 lmp_mpi -in lammps.in
```

### 2.2 Output Files

| File | Description |
|------|-------------|
| `dump.minimization` | Trajectory during minimization (every 400 steps) |
| `dump.Final` | Final relaxed structure (sorted by atom id) |
| `log.lammps` | LAMMPS log file |

### 2.3 Required LAMMPS Packages

Ensure your LAMMPS build includes:
- `MANYBODY` — for `sw/mod`, `tersoff`, `rebo`
- `INTERLAYER` — for `kolmogorov/crespi/z`, `ilp/graphene/hbn`

---

## Step 3: Layer Extraction

### 3.1 Create `cutpos.inp`

```
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
```

For orthorhombic TMD systems set `orthocell_12atom_sw = True`.

To cut to a smaller supercell (optional):
```
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
lattice_parameters = [3.28, 3.28, 35.0]
superlattice_vectors_block
4 0 0
0 4 0
0 0 1
```

### 3.2 Run cutpos

```bash
python ../../Src/cutpos.py
```

**Output:**
| File | Description |
|------|-------------|
| `relaxed_structure.cif` | Full relaxed structure |
| `layer_1.cif`, `layer_2.cif`, ... | Individual layer CIF files |
| `layer_1_coords.dat`, `layer_2_coords.dat`, ... | Coordinate files for analysis |

---

## Step 4: Analysis

### 4.1 Automated Analysis

```bash
python ../../Src/run_analysis.py
```

This sequentially runs:
1. `setup_analysis_dirs.py` — creates `InterlayerSpacingMap/` and `StrainMap/` directories and copies coordinate files
2. `generate_plot_inputs.py` — writes `input` config files for each analysis subdirectory
3. `plot_interlayer_spacing.py` — scatter and interpolated interlayer spacing maps (multi-layer only)
4. `plot_strains.py` — per-layer strain maps and histograms

### 4.2 Output Structure

```
MySystem/
├── InterlayerSpacingMap/        # (bilayer and multilayer only)
│   └── Layer_1-2/
│       ├── scatter_Layer_1-2.png
│       └── interpolated_Layer_1-2.png
└── StrainMap/
    ├── Layer_1/
    │   ├── strain_Layer_1.png
    │   └── hist_strain_Layer_1.png
    └── Layer_2/
        ├── strain_Layer_2.png
        └── hist_strain_Layer_2.png
```

> **Note for graphene**: `plot_strains.py` automatically detects single-species layers and filters to one z-sublayer before computing the strain KDTree, preventing spurious interlayer C-C distances from corrupting the strain calculation.

---

## Step 5: DFT Conversion (Optional)

### SIESTA Format

```bash
python ../../Src/cif2siesta.py relaxed_structure.cif output.fdf
```

### Quantum ESPRESSO Format

```bash
python ../../Src/cif2qe.py relaxed_structure.cif output.in
```

---

## Quick Reference

### Complete Workflow Commands

```bash
# 1. Generate structure
cd Examples/MySystem
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
| `potential file not found` | Copy potential files (`.KC`, `.rebo`, `.sw`, `.tersoff`) to working directory |
| `Tag mismatch` in cutpos | Verify tags in `layer*.inp` match what was used to generate `structure.lammps` |
| Wrong interlayer spacing | Check `translate_z` values in `layer*.inp` |
| `Neighbor list overflow` | Add `neigh_modify one 5000` after `neighbor` line in `lammps.in` |
| Strain map looks wrong (graphene) | Ensure `plot_strains.py` is up to date — single-species sublayer filter is applied automatically |
