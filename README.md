<p align="center">
  <img src="docs/images/logo.jpg" alt="TwisterASE Logo" width="400">
</p>

# TwisterASE

**A Python toolkit for generating and analyzing twisted layered material heterostructures**

TwisterASE builds supercells of twisted 2D materials (graphene, hBN, TMDs) and generates LAMMPS input files for molecular dynamics simulations with automatic interlayer potential configuration.

---

## Features

### Structure Generation
- ✅ **Twist angle control** via commensurate supercell indices (i-value or m,n values)
- ✅ **Multi-layer heterostructures** (graphene, hBN, TMDs)
- ✅ **Hexagonal and orthorhombic basis** support for TMDs
- ✅ **Automatic supercell construction** with periodic boundary conditions
- ✅ **CIF output** for visualization and DFT calculations

### LAMMPS Integration
- ✅ **Automatic LAMMPS input generation** (`structure.lammps`, `lammps.in`)
- ✅ **Material-specific potentials**:
  - Graphene: REBO, KC interlayer
  - hBN: Tersoff intralayer, KC interlayer
  - TMD: SW/mod intralayer, KC interlayer
  - Mixed systems: Separate potentials per layer
- ✅ **Interlayer interaction generation**:
  - Hexagonal TMD: Classification-based (metal-metal, chalcogen-chalcogen)
  - Orthorhombic TMD: Full 64 interactions per layer pair
  - hBN-TMD: All cross-layer pairs with KC potential
- ✅ **Atom style**: `atomic` (no charges/molecules)

### Post-Processing & Analysis
- ✅ **Layer extraction** from relaxed LAMMPS structures (`cutpos.py`)
- ✅ **Interlayer spacing analysis** with heatmaps and statistics
- ✅ **Strain analysis** with distribution plots
- ✅ **DFT converters**: CIF → SIESTA (.fdf), Quantum ESPRESSO (.in)

---

## Installation

### Requirements
```bash
# Core dependencies
pip install numpy scipy ase matplotlib
```

### Python Version
- Python 3.7+

### External Software (Optional)
- **LAMMPS** (for MD simulations)
- **SIESTA** or **Quantum ESPRESSO** (for DFT calculations)

---

## Quick Start

### 1. Generate Twisted Bilayer Structure

Create input files in your working directory:

**`twisterase.inp`** (main configuration):
```python
# Number of layers
n_layers = 2
hex_lattice = True

# Twist angle via commensurate index
i_value = 7  # ~21.8° twist angle

# Superlattice parameters
lattice_parameters = [2.46, 2.46, 35.0]
superlattice_vectors_block
1 0 0
0 1 0
0 0 1

# LAMMPS output
write_lammps = True
interlayer_potential = "C.KC"
```

**`layer1.inp`** (bottom layer):
```python
twist_angle = 0.0
lattice_parameters = [2.46, 2.46, 35.0]

lattice_vectors_block
1.0  0.0  0.0
-0.5  0.866025  0.0
0.0  0.0  1.0

start_atom_positions_block
C  0.0    0.0      0.0  1
C  0.333  0.666    0.0  2
end_atom_positions_block

intralayer_potential = "C.rebo"
```

**`layer2.inp`** (top layer):
```python
twist_angle = 21.787  # Will be overridden by i_value
lattice_parameters = [2.46, 2.46, 35.0]
translate_z = 3.35

lattice_vectors_block
1.0  0.0  0.0
-0.5  0.866025  0.0
0.0  0.0  1.0

start_atom_positions_block
C  0.0    0.0      0.0  3
C  0.333  0.666    0.0  4
end_atom_positions_block

intralayer_potential = "C.rebo"
```

### 2. Run Structure Generation

```bash
cd Test/MySystem
python ../../Src/twisterase.py
```

**Output:**
- `superlattice.cif` - Full twisted structure
- `structure.lammps` - LAMMPS data file
- `lammps.in` - LAMMPS input script

### 3. Run LAMMPS Simulation

```bash
lmp_serial -in lammps.in
```

**Note:** The structure file uses `atom_style full` with charges for ILP potential compatibility. Ensure your LAMMPS build includes:
- `pair_style sw/mod` (TMDs)
- `pair_style tersoff` (hBN)
- `pair_style kolmogorov/crespi/z` (interlayer)
- `pair_style rebo` (graphene)

---

## Supported Materials

### Graphene
- **Basis**: 2 atoms (hexagonal)
- **Intralayer**: REBO
- **Interlayer**: KC (Kolmogorov-Crespi)
- **Examples**: `Test/Graphene_Bilayer_Hex/`

### hBN (Hexagonal Boron Nitride)
- **Basis**: 2 atoms (B, N)
- **Intralayer**: Tersoff (separate instance per layer)
- **Interlayer**: KC
- **Examples**: `Test/hBN_Trilayer_Hex/`

### TMD (Transition Metal Dichalcogenides)
- **Hexagonal basis**: 3 atoms (1 metal + 2 chalcogens)
  - Intralayer: SW/mod
  - Interlayer: KC with classification (metal-metal, chalcogen-chalcogen)
  - Examples: `Test/WSe2_Bilayer_Hex/`
  
- **Orthorhombic basis**: 12 atoms (4 metals + 8 chalcogens)
  - Intralayer: SW/mod with numbered atoms (W1, W2, Se1, Se2, ...)
  - Interlayer: 64 interactions per layer pair
  - Examples: `Test/WS2_WSe2_Bilayer_Ortho/`

### Mixed Heterostructures
- **hBN + TMD**: Separate potentials, KC for all interlayer
  - Examples: `Test/tWSe2_hBN/`
- **Graphene + TMD**: Not yet implemented
- **hBN + Graphene**: Not yet implemented

---

## File Structure

```
TwisterASEv2/
├── Src/
│   ├── twisterase.py              # Main structure generator
│   ├── layer.py                   # Layer class
│   ├── file_io.py                 # Input file parsing
│   ├── transformations.py         # Twist angle calculations
│   ├── lammps_writers/
│   │   ├── structure_writer.py    # LAMMPS data file writer
│   │   ├── input_generator.py     # Main LAMMPS input generator
│   │   ├── input_generator_hbn_tmd.py  # Mixed hBN+TMD generator
│   │   └── material_detector.py   # Material type detection
│   ├── cutpos.py                  # Layer extraction tool
│   ├── analysis.py                # Structure analysis
│   ├── cif2siesta.py              # CIF → SIESTA converter
│   ├── cif2qe.py                  # CIF → Quantum ESPRESSO converter
│   └── run_analysis.py            # Automated analysis workflow
│
├── Test/
│   ├── Graphene_Bilayer_Hex/      # Twisted bilayer graphene
│   ├── hBN_Trilayer_Hex/          # hBN trilayer
│   ├── WSe2_Bilayer_Hex/          # Hexagonal WSe2 bilayer
│   ├── WS2_WSe2_Bilayer_Ortho/    # Orthorhombic WS2/WSe2
│   └── tWSe2_hBN/                 # Mixed TMD-hBN heterostructure
│
├── docs/
│   ├── images/                    # Logo and other images
│   ├── input_files.md             # Input file format reference
│   ├── workflow.md                # Complete workflow guide
│   ├── cutpos.md                  # Layer extraction documentation
│   └── converters.md              # DFT converter documentation
│
└── README.md                      # This file
```

---

## Workflow

### Complete Workflow: Structure → Simulation → Analysis

```bash
# 1. Generate structure
cd Test/MySystem
python ../../Src/twisterase.py

# 2. Run LAMMPS relaxation
lmp_serial -in lammps.in

# 3. Extract individual layers
python ../../Src/cutpos.py

# 4. Run analysis
python ../../Src/run_analysis.py

# 5. Convert to DFT format (optional)
python ../../Src/cif2siesta.py layer_1.cif layer_1.fdf
python ../../Src/cif2qe.py layer_1.cif layer_1.in
```

---

## Input File Reference

### `twisterase.inp` Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_layers` | int | Number of layers |
| `hex_lattice` | bool | Use hexagonal lattice (True) or custom |
| `i_value` | int | Commensurate supercell index (hex only) |
| `mn_values` | [int, int] | Alternative twist angle specification |
| `lattice_parameters` | [float, float, float] | Lattice constants [a, b, c] |
| `superlattice_vectors_block` | 3×3 matrix | Superlattice vectors |
| `write_lammps` | bool | Generate LAMMPS files |
| `interlayer_potential` | str | Interlayer potential file (e.g., "C.KC") |

### `layer*.inp` Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `twist_angle` | float | Rotation angle (degrees) |
| `lattice_parameters` | [float, float, float] | Layer lattice constants |
| `translate_z` | float | Vertical offset (Å) |
| `lattice_vectors_block` | 3×3 matrix | Unit cell vectors |
| `start_atom_positions_block` | list | Atom positions (symbol x y z tag) |
| `intralayer_potential` | str | Potential file (e.g., "tmd.sw") |
| `orthocell_transformation_matrix` | 3×3 matrix | Ortho basis transform (TMD only) |
| `start_atom_positions_ortho_block` | list | Ortho basis atoms (TMD only) |

---

## LAMMPS Output Format

### For Mixed hBN+TMD Systems

The generator produces `lammps.in` with:

```lammps
# General
units           metal
dimension       3
atom_style      atomic
neighbor        0.3 bin

# Structure
boundary        p p p
read_data       structure.lammps
mass 1 183.84   # W
mass 2 78.971   # Se
...

# potential definitions
pair_style hybrid/overlay sw/mod sw/mod tersoff tersoff kolmogorov/crespi/z 14.0 ...
# Number of interlayer interactions = 8

# intralayer interactions
pair_coeff * * sw/mod 1 tmd.sw W Se Se NULL NULL NULL NULL NULL NULL NULL
pair_coeff * * sw/mod 2 tmd.sw NULL NULL NULL W Se Se NULL NULL NULL NULL
pair_coeff * * tersoff 1 BNC.tersoff NULL NULL NULL NULL NULL NULL B N NULL NULL
pair_coeff * * tersoff 2 BNC.tersoff NULL NULL NULL NULL NULL NULL NULL NULL B N

pair_coeff * * lj/cut 0.0 3.4

# interlayer interaction
pair_coeff 7 2 kolmogorov/crespi/z 1 SeBN.KC NULL Se NULL NULL NULL NULL B NULL NULL NULL
...

# Optimize at 0 K
dump            1 all custom 400 dump.minimization id type x y z
thermo          500
min_style       fire
minimize        0.0 1.0e-4 1000000 1000000
write_dump all atom dump.Final
print "Done!"
```

---

## Known Issues & Limitations

### Current Limitations
1. **Structure file format**: Uses `atom_style full` with charges, but `lammps.in` uses `atom_style atomic`
   - **Impact**: LAMMPS fails to read structure file
   - **Status**: Needs structure writer update to match atomic style
   
2. **Mixed systems**: Only hBN+TMD fully implemented
   - Graphene+TMD: Not implemented
   - hBN+Graphene: Not implemented

3. **Layer ordering**: Structure file processes layers sequentially
   - hBN layers get remapped types based on position
   - May not match expected ordering in some cases

### Workarounds
- For structure file format issue: Manually edit `structure.lammps` to remove molecule IDs and charges
- For layer ordering: Rename layer files to control processing order

---

## Examples

See `Test/` directory for working examples:

- **`Graphene_Bilayer_Hex/`**: Twisted bilayer graphene (21.8°)
- **`hBN_Trilayer_Hex/`**: Three-layer hBN stack
- **`WSe2_Bilayer_Hex/`**: Hexagonal WSe2 bilayer
- **`WS2_WSe2_Bilayer_Ortho/`**: Orthorhombic WS2/WSe2 heterostructure
- **`tWSe2_hBN/`**: Mixed TMD-hBN system (4 layers)

Each example includes:
- Input files (`twisterase.inp`, `layer*.inp`, `cutpos.inp`)
- Potential files (`.KC`, `.sw`, `.tersoff`, `.rebo`)
- Generated structures (`.cif`, `.lammps`)

---

## Contributing

TwisterASEv2 is under active development. Current priorities:

1. Fix structure file format to use `atom_style atomic`
2. Implement Graphene+TMD mixed systems
3. Add support for more TMD materials (MoS2, MoSe2, WS2)
4. Improve documentation and examples

---

## References

### Potentials
- **KC (Kolmogorov-Crespi)**: Interlayer van der Waals interactions
- **REBO**: Reactive empirical bond order (graphene)
- **Tersoff**: Three-body potential (hBN)
- **SW/mod**: Modified Stillinger-Weber (TMDs)

### Related Tools
- **ASE**: Atomic Simulation Environment (structure manipulation)
- **LAMMPS**: Large-scale Atomic/Molecular Massively Parallel Simulator

---

## License

[Add license information]

---

## Contact

[Add contact information]
