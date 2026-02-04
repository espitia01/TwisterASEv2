# DFT Converters

Convert CIF structures to DFT input formats for SIESTA and Quantum ESPRESSO.

---

## CIF to SIESTA

### Usage

```bash
python ../../Src/cif2siesta.py <input.cif> <output.fdf>
```

### Example

```bash
python ../../Src/cif2siesta.py layer_1.cif layer_1.fdf
```

### Output Format

The converter generates a SIESTA `.fdf` file with:
- Lattice vectors
- Atomic coordinates (fractional)
- Species labels and atomic numbers

---

## CIF to Quantum ESPRESSO

### Usage

```bash
python ../../Src/cif2qe.py <input.cif> <output.in>
```

### Example

```bash
python ../../Src/cif2qe.py layer_1.cif layer_1.in
```

### Output Format

The converter generates a QE input file with:
- `&SYSTEM` namelist with cell parameters
- `CELL_PARAMETERS` block
- `ATOMIC_POSITIONS` block (crystal coordinates)

---

## Supported Elements

Both converters support all elements present in typical 2D materials:
- **Graphene**: C
- **hBN**: B, N
- **TMDs**: Mo, W, S, Se, Te

---

## Notes

- Converters read the full unit cell from the CIF file
- Atomic positions are converted to fractional coordinates
- Additional DFT parameters (k-points, cutoffs, etc.) must be added manually
