# Layer Extraction Tool (cutpos.py)

Extract individual layers from relaxed LAMMPS structures for analysis.

---

## Overview

`cutpos.py` reads a LAMMPS dump file and extracts individual layers based on atom type tags defined in the layer input files.

---

## Requirements

Before running `cutpos.py`, ensure you have:

1. **`cutpos.inp`** - Configuration file
2. **`layer*.inp`** - Layer definition files (same as used for structure generation)
3. **`dump.Final`** - LAMMPS dump file with relaxed coordinates

---

## Configuration: cutpos.inp

```python
n_layers = 2
lammps_dump = "dump.Final"
orthocell_12atom_sw = False
# cut = 15.0  # Optional: z-coordinate to cut structure
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_layers` | int | required | Number of layers to extract |
| `lammps_dump` | str | "dump.Final" | LAMMPS dump file to read |
| `orthocell_12atom_sw` | bool | False | Use 12-atom orthorhombic basis |
| `cut` | float | None | Optional z-coordinate cutoff |

---

## Usage

```bash
python ../../Src/cutpos.py
```

---

## Output Files

| File | Description |
|------|-------------|
| `relaxed_structure.cif` | Full relaxed structure |
| `layer_N.cif` | Individual layer N as CIF |
| `layer_N_coords.dat` | Layer N coordinates (x, y, z format) |

---

## How It Works

1. **Load Configuration**: Reads `cutpos.inp` and `layer*.inp` files
2. **Build Tag Mapping**: Maps atom type tags to layers based on layer definitions
3. **Read Dump File**: Parses LAMMPS dump file for atom positions and types
4. **Extract Layers**: Separates atoms by type into individual layers
5. **Write Output**: Generates CIF and coordinate files for each layer

---

## Tag Mapping

The tool uses tags from layer input files to identify which atoms belong to which layer:

```
Layer 1: tags [1, 2, 3] -> atoms with types 1, 2, 3
Layer 2: tags [4, 5, 6] -> atoms with types 4, 5, 6
```

**Important**: Ensure tags in layer files match the atom types in your LAMMPS structure.

---

## Example Workflow

```bash
# 1. Ensure you have the required files
ls
# cutpos.inp  layer1.inp  layer2.inp  dump.Final  ...

# 2. Run extraction
python ../../Src/cutpos.py

# 3. Check output
ls *.cif *.dat
# relaxed_structure.cif  layer_1.cif  layer_2.cif
# layer_1_coords.dat  layer_2_coords.dat
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Tag validation failed" | Check that tags in layer files match structure atom types |
| Missing atoms in layer | Verify tag assignments are unique and complete |
| Wrong layer assignment | Ensure layer files are numbered correctly (layer1, layer2, ...) |
