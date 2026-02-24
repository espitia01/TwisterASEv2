# Input File Reference

This document describes the format and parameters for all TwisterASE input files.

---

## Main Configuration: `twisterase.inp`

Controls overall structure generation. All scripts (`twisterase.py`, `cutpos.py`, `run_analysis.py`) must be run from the directory containing this file.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_layers` | int | Number of layers in the heterostructure |
| `hex_lattice` | bool | `True` for hexagonal geometry (required for `i_value`/`mn_values`) |
| `i_value` | int | Commensurate supercell index — sets twist angle for layers 1 & 2 |
| `mn_values` | [int, int] | Alternative twist angle via [m, n] integers |
| `lattice_parameters` | [float, float, float] | Lattice constants [a, b, c] in Å — used with `superlattice_vectors_block` |
| `superlattice_vectors_block` | 3×3 | Superlattice vectors in units of `lattice_parameters` |
| `write_lammps` | bool | Generate `structure.lammps` and `lammps.in` |
| `interlayer_potential` | str | KC potential filename (e.g. `"CC.KC"`, `"MoWSSe.KC"`) |

> **Twist angle**: Use either `i_value` or `mn_values` (not both). When set, the superlattice vectors are computed automatically from the layer 1 unit cell and `lattice_parameters`/`superlattice_vectors_block` are not needed. Only layers 1 and 2 get their twist angles overridden; layers 3+ keep the `twist_angle` from their own file.

### Example — Twisted Bilayer Graphene (`i_value`)

```
n_layers = 2
hex_lattice = True
i_value = 9

write_lammps = True
interlayer_potential = "CC.KC"
```

### Example — Custom Superlattice Vectors

```
n_layers = 2
hex_lattice = True

lattice_parameters = [3.3, 3.3, 35.0]
superlattice_vectors_block
10.0  0.0  0.0
0.0  10.0  0.0
0.0   0.0  1.0

write_lammps = True
interlayer_potential = "MoWSSe.KC"
```

---

## Layer Configuration: `layer*.inp`

Each layer requires a separate input file (`layer1.inp`, `layer2.inp`, etc.).

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `twist_angle` | float | Rotation angle in degrees (overridden by `i_value`/`mn_values` for layers 1 & 2) |
| `lattice_parameters` | [float, float, float] | Layer lattice constants [a, b, c] in Å |
| `translate_z` | float | Vertical offset in Å — sets interlayer spacing |
| `lattice_vectors_block` | 3×3 | Unit cell vectors in fractional coordinates (rows = a1, a2, a3) |
| `start_atom_positions_block` | block | Atom positions in scaled coordinates |
| `end_atom_positions_block` | — | Closes the atom block |
| `intralayer_potential` | str | Potential filename for this layer |
| `orthocell_transformation_matrix` | 3×3 | Transform from hex to orthorhombic basis (TMD ortho only) |
| `start_atom_positions_ortho_block` | block | Orthorhombic basis atom positions (TMD ortho only) |
| `end_atom_positions_ortho_block` | — | Closes the ortho atom block |

### Atom Position Format

```
start_atom_positions_block
Symbol  x  y  z  tag
...
end_atom_positions_block
```

- **Symbol**: Element symbol (e.g. `C`, `W`, `Se`, `B`, `N`)
- **x, y, z**: Fractional (scaled) coordinates
- **tag**: Unique integer atom type identifier across all layers

### Example — Graphene Layer

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

### Example — TMD Layer (Hexagonal, MoSe2)

```
twist_angle = 0.0
lattice_parameters = [3.3, 3.3, 35.0]

lattice_vectors_block
1.0 0.0 0.0
0.5 0.8660254038 0.0
0.0 0.0 1.0

start_atom_positions_block
Mo  0.0          0.0          0.168272109  1
Se  0.666666667  0.666666667  0.124125946  2
Se  0.666666667  0.666666667  0.212418272  3
end_atom_positions_block

translate_z = 0.0
intralayer_potential = "tmd.sw"
```

### Example — hBN Layer

```
twist_angle = 0.0
lattice_parameters = [2.504, 2.504, 35.0]

lattice_vectors_block
1.0 0.0 0.0
0.5 0.8660254038 0.0
0.0 0.0 1.0

start_atom_positions_block
B  0.0          0.0          0.026  7
N  0.666666667  0.666666667  0.026  8
end_atom_positions_block

translate_z = 0.0
intralayer_potential = "BNC.tersoff"
```

---

## Post-Processing: `cutpos.inp`

Configuration for layer extraction from relaxed LAMMPS dump files.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_layers` | int | required | Number of layers to extract |
| `lammps_dump` | str | `"dump.Final"` | LAMMPS dump file to process |
| `orthocell_12atom_sw` | bool | `False` | Use 12-atom orthorhombic basis for layer assignment |
| `lattice_parameters` | [float, float, float] | None | Lattice constants for optional cut cell |
| `superlattice_vectors_block` | 3×3 | None | Cut to smaller supercell (optional) |

### Example — Basic

```
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
```

### Example — With Supercell Cut

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

---

## Atom Type Tags

Tags identify atom types across layers and must be **unique** for each distinct atom type in the system. They are used to:
- Assign LAMMPS atom types in `structure.lammps`
- Map atoms back to layers in `cutpos.py`
- Generate correct `pair_coeff` lines in `lammps.in`

### Tag Conventions by System

**Graphene bilayer** (2 atoms/layer):
- Layer 1: C→1, C→2
- Layer 2: C→3, C→4

**TMD bilayer, hexagonal** (3 atoms/layer):
- Layer 1: Mo→1, Se→2, Se→3
- Layer 2: Mo→4, Se→5, Se→6

**TMD bilayer, orthorhombic** (12 atoms/layer):
- Layer 1: tags 1–12
- Layer 2: tags 13–24

**hBN + TMD (hex), 4-layer stack**:
- TMD Layer 1: tags 1–3; TMD Layer 2: tags 4–6
- hBN below: B→7, N→8; hBN above: B→9, N→10

**hBN + TMD (ortho), 4-layer stack**:
- TMD Layer 1: tags 1–12; TMD Layer 2: tags 13–24
- hBN below: B→25, N→26; hBN above: B→27, N→28
