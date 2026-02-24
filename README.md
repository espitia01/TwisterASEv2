<p align="center">
  <img src="docs/images/logo.png" alt="TwisterASE Logo">
</p>

# TwisterASE

**A Python toolkit for generating and analyzing twisted layered material heterostructures**

TwisterASE builds supercells of twisted 2D materials (graphene, hBN, TMDs) and generates LAMMPS input files for molecular dynamics simulations with automatic interlayer potential configuration.

---

## Features

### Structure Generation
- **Twist angle control** via commensurate supercell indices (i-value or m,n values)
- **Multi-layer heterostructures** (graphene, hBN, TMDs)
- **Hexagonal and orthorhombic basis** support for TMDs
- **Automatic supercell construction** with periodic boundary conditions
- **CIF output** for visualization and DFT calculations

### LAMMPS Integration
- **Automatic LAMMPS input generation** (`structure.lammps`, `lammps.in`)
- **Material-specific potentials**:
  - Graphene: REBO, KC interlayer
  - hBN: Tersoff intralayer, KC interlayer
  - TMD: SW/mod intralayer, KC interlayer
  - Mixed systems: Separate potentials per layer
- **Interlayer interaction generation**:
  - Hexagonal TMD: Classification-based (metal-metal, chalcogen-chalcogen)
  - Orthorhombic TMD: Full 64 interactions per layer pair
  - hBN-TMD: All cross-layer pairs with KC potential
- **Atom style**: `atomic` (no charges/molecules)

### Post-Processing & Analysis
- **Layer extraction** from relaxed LAMMPS structures (`cutpos.py`)
- **Interlayer spacing analysis** with heatmaps and statistics
- **Strain analysis** with distribution plots
- **DFT converters**: CIF → SIESTA (.fdf), Quantum ESPRESSO (.in)

---

## Installation

### Conda (Recommended)

```bash
conda env create -f environment.yml
conda activate twisterase
```

### pip

```bash
pip install numpy scipy ase matplotlib
```

### Python Version
- Python 3.9+

### External Software (Optional)
- **LAMMPS** (for MD simulations)
- **SIESTA** or **Quantum ESPRESSO** (for DFT calculations)

---

## Worked Example: Twisted Bilayer Graphene

This reproduces `Examples/Graphene_Bilayer_Hex/` — a commensurate twisted bilayer graphene at `i_value = 9` (~13.17°).

### Step 1 — Input files

**`twisterase.inp`**
```
n_layers = 2
hex_lattice = True
i_value = 9

write_lammps = True
interlayer_potential = "CC.KC"
```

**`layer1.inp`** — bottom layer (tags 1, 2)
```
twist_angle = 0.0
lattice_parameters = [2.46, 2.46, 35.0]

lattice_vectors_block
1.0 0.0 0.0
0.5 0.8660254038 0.0
0.0 0.0 1.0

start_atom_positions_block
C   0.000000000   0.000000000   0.5  1
C   0.666666666   0.666666666   0.5  2
end_atom_positions_block

translate_z = 0.0
intralayer_potential = "CH.rebo"
```

**`layer2.inp`** — top layer (tags 3, 4; twist set automatically by `i_value`)
```
twist_angle = 0.0
lattice_parameters = [2.46, 2.46, 35.0]

lattice_vectors_block
1.0 0.0 0.0
0.5 0.8660254038 0.0
0.0 0.0 1.0

start_atom_positions_block
C   0.000000000   0.000000000   0.5  3
C   0.666666666   0.666666666   0.5  4
end_atom_positions_block

translate_z = 3.3
intralayer_potential = "CH.rebo"
```

> **Tag rule**: each atom type must have a unique integer tag across all layers.

### Step 2 — Generate structure

```bash
cd Examples/Graphene_Bilayer_Hex
python ../../Src/twisterase.py
```

Produces: `superlattice.cif`, `structure.lammps`, `lammps.in`

### Step 3 — Generated `lammps.in`

```lammps
units           metal
dimension       3
atom_style      atomic
neighbor        0.3 bin

boundary        p p p
read_data       structure.lammps
mass 1 12.011
mass 2 12.011
mass 3 12.011
mass 4 12.011

pair_style      hybrid/overlay rebo kolmogorov/crespi/z 14.0

pair_coeff      * * rebo CH.rebo C C C C
pair_coeff      1 3 kolmogorov/crespi/z CC.KC C C C C

dump            1 all custom 400 dump.minimization id type x y z
min_style       fire
minimize        0.0 1.0e-8 1000000 1000000
undump          1
write_dump      all custom dump.Final id type x y z modify sort id
print "Done!"
```

### Step 4 — Run LAMMPS

```bash
lmp_serial -in lammps.in
# or in parallel:
mpirun -np 8 lmp_mpi -in lammps.in
```

### Step 5 — Extract layers

**`cutpos.inp`**:
```
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
```

```bash
python ../../Src/cutpos.py
```

Produces: `relaxed_structure.cif`, `layer_1.cif`, `layer_2.cif`, `layer_1_coords.dat`, `layer_2_coords.dat`

### Step 6 — Run analysis

```bash
python ../../Src/run_analysis.py
```

Produces interlayer spacing maps and per-layer strain maps in `InterlayerSpacingMap/` and `StrainMap/`.

---

## Supported Materials

| Material | Basis | Intralayer | Interlayer | Examples |
|----------|-------|-----------|-----------|---------|
| **Graphene** | 2 C atoms (hex) | REBO (`CH.rebo`) | KC (`CC.KC`) — one `pair_coeff` per adjacent pair | `Graphene_Bilayer_Hex/` |
| **hBN** | 2 atoms B, N (hex) | Tersoff (`BNC.tersoff`) per layer | ILP (`BNCH.ILP`) | `Bilayer_hBN/`, `hBN_Trilayer_Hex/` |
| **TMD hex** | 3 atoms: TM + 2 X (hex) | SW/mod (`tmd.sw`) | KC — z-classified (TM, X_l, X_u) | `tMoSe2/`, `WSe2_Bilayer_Hex/` |
| **TMD ortho** | 12 atoms: 4 TM + 8 X (ortho) | SW/mod | KC — 64 interactions per layer pair | `tMoSe2_Ortho/`, `WS2_WSe2_Bilayer_Ortho/` |
| **hBN + TMD** | mixed | SW/mod + Tersoff per layer | KC chalcogen-specific (`SeBN.KC`, `SBN.KC`) | `tWSe2_hBN/`, `tMoSe2_Ortho_hBN/` |

> Graphene+TMD and hBN+Graphene mixed systems are not yet implemented.

---

## File Structure

```
TwisterASEv2/
├── Src/
│   ├── twisterase.py              # Main structure generator
│   ├── layer.py                   # Layer class
│   ├── file_io.py                 # Input file parsing
│   ├── transformations.py         # Twist angle & geometry utilities
│   ├── lammps_writers/
│   │   ├── structure_writer.py    # LAMMPS data file writer
│   │   ├── input_generator.py     # LAMMPS input generator (graphene, hBN, TMD)
│   │   ├── input_generator_hbn_tmd.py  # Mixed hBN+TMD generator
│   │   └── material_detector.py   # Material type detection
│   ├── cutpos.py                  # Layer extraction from dump.Final
│   ├── run_analysis.py            # Automated analysis workflow
│   ├── setup_analysis_dirs.py     # Analysis directory setup
│   ├── generate_plot_inputs.py    # Plot input file generator
│   ├── plot_interlayer_spacing.py # Interlayer spacing heatmaps
│   ├── plot_strains.py            # Strain maps and histograms
│   ├── cif2siesta.py              # CIF → SIESTA converter
│   └── cif2qe.py                  # CIF → Quantum ESPRESSO converter
│
├── Examples/
│   ├── Graphene_Bilayer_Hex/      # Twisted bilayer graphene (i_value=9, ~13.2°)
│   ├── Bilayer_hBN/               # hBN bilayer
│   ├── hBN_Trilayer_Hex/          # hBN trilayer
│   ├── WSe2_Bilayer_Hex/          # Hexagonal WSe2 bilayer
│   ├── WS2_WSe2_Bilayer_Ortho/    # Orthorhombic WS2/WSe2 heterostructure
│   ├── tMoSe2/                    # Twisted MoSe2 (hex, i_value=5)
│   ├── tMoSe2_Ortho/              # Twisted MoSe2 (orthorhombic basis)
│   ├── tMoSe2_Ortho_hBN/          # MoSe2 ortho + hBN encapsulation
│   ├── tWSe2_hBN/                 # Mixed TMD-hBN heterostructure
│   └── MoSe2-10x10_SC/            # Large MoSe2 supercell
│
├── docs/
│   ├── input_files.md             # Full input file format reference
│   ├── workflow.md                # Complete workflow guide
│   └── images/
│
├── environment.yml
└── README.md
```

---

## Input File Reference

### `twisterase.inp`

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_layers` | int | Number of layers |
| `hex_lattice` | bool | `True` for hexagonal geometry |
| `i_value` | int | Commensurate index — sets twist angle for layers 1 & 2 automatically |
| `mn_values` | [int, int] | Alternative: twist angle via [m, n] |
| `lattice_parameters` | [float, float, float] | Lattice constants [a, b, c] in Å (used with `superlattice_vectors_block`) |
| `superlattice_vectors_block` | 3×3 | Superlattice vectors in units of `lattice_parameters` |
| `write_lammps` | bool | Generate `structure.lammps` and `lammps.in` |
| `interlayer_potential` | str | KC potential filename (e.g. `"CC.KC"`, `"MoWSSe.KC"`) |

> When `i_value` or `mn_values` is set, superlattice vectors are computed automatically from the layer 1 unit cell. Only layers 1 and 2 have their twist angles overridden; additional layers use the `twist_angle` from their own file.

### `layer*.inp`

| Parameter | Type | Description |
|-----------|------|-------------|
| `twist_angle` | float | Rotation angle in degrees |
| `lattice_parameters` | [float, float, float] | Layer lattice constants [a, b, c] in Å |
| `translate_z` | float | Vertical offset in Å (sets interlayer spacing) |
| `lattice_vectors_block` | 3×3 | Unit cell vectors in fractional coordinates |
| `start_atom_positions_block` | block | `Symbol x y z tag` in scaled coordinates |
| `intralayer_potential` | str | Potential filename for this layer |
| `orthocell_transformation_matrix` | 3×3 | Hex → ortho transform (TMD ortho only) |
| `start_atom_positions_ortho_block` | block | Orthorhombic basis atom positions (TMD ortho only) |

### `cutpos.inp`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_layers` | int | required | Number of layers to extract |
| `lammps_dump` | str | `"dump.Final"` | LAMMPS dump file |
| `orthocell_12atom_sw` | bool | `False` | Use orthorhombic basis for layer assignment |
| `lattice_parameters` | [float, float, float] | None | Lattice constants for optional cut cell |
| `superlattice_vectors_block` | 3×3 | None | Cut to smaller supercell (optional) |

---

## Atom Type Tags

Each atom type must have a **unique integer tag** across all layers. Tags map LAMMPS atom types to chemical species and layers.

**Graphene bilayer** (tags 1–4):
- Layer 1: C→1, C→2
- Layer 2: C→3, C→4

**TMD bilayer, hexagonal** (tags 1–6):
- Layer 1: Mo→1, Se→2, Se→3
- Layer 2: Mo→4, Se→5, Se→6

**TMD bilayer, orthorhombic** (tags 1–24):
- Layer 1: 12 atoms, tags 1–12
- Layer 2: 12 atoms, tags 13–24

**hBN + TMD** (hBN tags follow max TMD tag):
- TMD Layer 1 (hex): tags 1–3; TMD Layer 2: tags 4–6
- hBN below: B→7, N→8; hBN above: B→9, N→10

---

## Known Limitations

1. **Mixed systems**: Only hBN+TMD is fully implemented. Graphene+TMD and hBN+Graphene are not yet supported.
2. **Twist angle override**: `i_value`/`mn_values` only sets the twist for layers 1 and 2. Layers 3+ use `twist_angle` from their input file.
3. **Ortho z-tolerance**: Chalcogen partitioning uses a fixed z-tolerance. Highly corrugated post-relaxation structures may need manual verification.

---

## Complete Workflow (Quick Reference)

```bash
# 1. Generate structure and LAMMPS files
cd Examples/MySystem
python ../../Src/twisterase.py

# 2. Run LAMMPS relaxation
lmp_serial -in lammps.in
# or: mpirun -np 8 lmp_mpi -in lammps.in

# 3. Extract individual layers from relaxed dump
python ../../Src/cutpos.py

# 4. Run analysis (interlayer spacing + strain maps)
python ../../Src/run_analysis.py

# 5. Convert to DFT format (optional)
python ../../Src/cif2siesta.py relaxed_structure.cif output.fdf
python ../../Src/cif2qe.py relaxed_structure.cif output.in
```

See `docs/workflow.md` for the full step-by-step guide and `docs/input_files.md` for complete parameter reference.
