# Input File Reference

This document describes the format and parameters for TwisterASE input files.

---

## Main Configuration: `twisterase.inp`

The main configuration file controls the overall structure generation.

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_layers` | int | Number of layers in the heterostructure |

### Twist Angle Specification

Choose **one** of the following methods:

| Parameter | Type | Description |
|-----------|------|-------------|
| `i_value` | int | Commensurate supercell index (hexagonal lattice only) |
| `mn_values` | [int, int] | Alternative twist angle specification [m, n] |

### Lattice Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `hex_lattice` | bool | Use hexagonal lattice (True) or custom |
| `lattice_parameters` | [float, float, float] | Lattice constants [a, b, c] in Å |
| `superlattice_vectors_block` | 3×3 matrix | Superlattice transformation vectors |

### Output Options

| Parameter | Type | Description |
|-----------|------|-------------|
| `write_lammps` | bool | Generate LAMMPS structure and input files |
| `interlayer_potential` | str | Interlayer potential file (e.g., "C.KC", "MoWSSe.KC") |

### Example

```python
n_layers = 2
hex_lattice = True

i_value = 7  # ~21.8° twist angle

lattice_parameters = [3.28, 3.28, 35.0]

superlattice_vectors_block
1 0 0
0 1 0
0 0 1

write_lammps = True
interlayer_potential = "WSe2.KC"
```

---

## Layer Configuration: `layer*.inp`

Each layer requires a separate input file (`layer1.inp`, `layer2.inp`, etc.).

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `twist_angle` | float | Rotation angle in degrees (may be overridden by i_value) |
| `lattice_parameters` | [float, float, float] | Layer lattice constants [a, b, c] |
| `lattice_vectors_block` | 3×3 matrix | Unit cell vectors in fractional coordinates |
| `start_atom_positions_block` | list | Atom positions with optional tags |

### Optional Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `translate_z` | float | Vertical offset in Å |
| `strain_percentage` | [float, float, float] | Applied strain [εx, εy, εz] |
| `intralayer_potential` | str | Potential file for this layer |

### Orthorhombic Basis (TMD only)

| Parameter | Type | Description |
|-----------|------|-------------|
| `orthocell_transformation_matrix` | 3×3 matrix | Transformation to orthorhombic basis |
| `start_atom_positions_ortho_block` | list | Orthorhombic basis atom positions |

### Atom Position Format

```
start_atom_positions_block
Symbol  x  y  z  [tag]
end_atom_positions_block
```

- **Symbol**: Element symbol (e.g., W, Se, B, N, C)
- **x, y, z**: Fractional coordinates
- **tag** (optional): Integer atom type identifier

### Examples

#### Graphene Layer
```python
twist_angle = 0.0
lattice_parameters = [2.46, 2.46, 35.0]

lattice_vectors_block
1.0   0.0       0.0
-0.5  0.866025  0.0
0.0   0.0       1.0

start_atom_positions_block
C  0.0    0.0    0.0  1
C  0.333  0.666  0.0  2
end_atom_positions_block

intralayer_potential = "C.rebo"
```

#### TMD Layer (Hexagonal)
```python
twist_angle = 0.0
lattice_parameters = [3.28, 3.28, 35.0]
translate_z = 0.0

lattice_vectors_block
-0.01990469  0.99980188  0.0
-0.87580618  0.48266298  0.0
0.0  0.0  1.0

start_atom_positions_block
W   0.0          0.0          0.04206802725  1
Se  0.333333333  0.333333333  0.02546802725  2
Se  0.333333333  0.333333333  0.05866802725  3
end_atom_positions_block

intralayer_potential = "tmd.sw"
```

#### hBN Layer
```python
twist_angle = 0.0
lattice_parameters = [2.4967, 2.4967, 100.0]

lattice_vectors_block
-0.01990469  0.99980188  0.0
-0.87580618  0.48266298  0.0
0.0  0.0  1.0

start_atom_positions_block
B  0.0          0.0          0.026  7
N  0.666666667  0.666666667  0.026  8
end_atom_positions_block

intralayer_potential = "BNC.tersoff"
```

---

## Post-Processing: `cutpos.inp`

Configuration for layer extraction from relaxed structures.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_layers` | int | Number of layers to extract |
| `lammps_dump` | str | LAMMPS dump file to process (default: "dump.Final") |
| `orthocell_12atom_sw` | bool | Use 12-atom orthorhombic basis |
| `cut` | float | Optional z-coordinate to cut structure |

### Example

```python
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
```

---

## Atom Type Tags

### Importance of Tags

Tags identify atom types across layers and must be **unique** for each distinct atom type in the system.

### Convention

For a 4-layer hBN-TMD-TMD-hBN system:
- **TMD Layer 1**: tags 1, 2, 3 (W, Se_lower, Se_upper)
- **TMD Layer 2**: tags 4, 5, 6 (W, Se_lower, Se_upper)
- **hBN Layer 1**: tags 7, 8 (B, N) - auto-assigned
- **hBN Layer 2**: tags 9, 10 (B, N) - auto-assigned

**Note**: hBN layer tags are automatically assigned after all TMD types.
