"""
cutpos.py - Extract individual layers from LAMMPS dump.Final

Reads the relaxed structure from dump.Final, cuts it to a smaller supercell,
and extracts individual layer positions as CIF and coordinate files.
"""

import numpy as np
import sys
import os
from ase.io import read, write
from ase import Atoms
from layer import Layer
from transformations import hexcut


def parse_cutpos_input(filename='cutpos.inp'):
    """
    Parse cutpos.inp file.
    
    Expected format:
        n_layers = 2
        lammps_dump = dump.Final
        orthocell_12atom_sw = False
        lattice_parameters = [3.16, 3.16, 35.0]
        superlattice_vectors_block
        1 0 0
        0 1 0
        0 0 1
    """
    config = {
        'n_layers': None,
        'lammps_dump': 'dump.Final',
        'orthocell_12atom_sw': False,
        'lattice_parameters': None,
        'cut_vectors': None
    }
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"[FAIL] {filename} not found! Using defaults.")
        return config
    
    for iline, line in enumerate(lines):
        keyword, _, value = line.partition('#')[0].partition('=')
        keyword = keyword.strip()
        
        if "n_layers" in keyword:
            config['n_layers'] = eval(value)
        elif "lammps_dump" in keyword:
            config['lammps_dump'] = eval(value)
        elif "orthocell_12atom_sw" in keyword:
            config['orthocell_12atom_sw'] = eval(value)
        elif "lattice_parameters" in keyword:
            config['lattice_parameters'] = np.array(eval(value))
        elif "superlattice_vectors_block" in keyword:
            cut_vecs = np.zeros((3, 3))
            for iv in range(3):
                cut_vecs[iv] = [float(x) for x in lines[iline+1+iv].split()]
            config['cut_vectors'] = cut_vecs
    
    return config


def parse_lammps_dump(dump_file):
    """
    Parse LAMMPS dump file to extract atom types.
    Returns: atoms (ASE Atoms), lammps_types (array of LAMMPS types)
    """
    with open(dump_file, 'r') as f:
        lines = f.readlines()
    
    # Find ITEM: ATOMS line
    atoms_line_idx = None
    for i, line in enumerate(lines):
        if 'ITEM: ATOMS' in line:
            atoms_line_idx = i
            break
    
    if atoms_line_idx is None:
        raise ValueError("ITEM: ATOMS not found in dump file")
    
    # Parse atom data: id type x y z (or id type xs ys zs)
    atom_data = []
    for i in range(atoms_line_idx + 1, len(lines)):
        line = lines[i].strip()
        if not line or 'ITEM:' in line:
            break
        parts = line.split()
        if len(parts) >= 5:
            atom_id = int(parts[0])
            atom_type = int(parts[1])
            atom_data.append((atom_id, atom_type))
    
    # Sort by atom id to match ASE ordering
    atom_data.sort(key=lambda x: x[0])
    lammps_types = np.array([atype for _, atype in atom_data])
    
    # Read structure with ASE
    atoms = read(dump_file, format='lammps-dump-text', index=-1)
    
    return atoms, lammps_types


def main():
    print(f"\n{'='*60}")
    print("CUTPOS: Extract Layers from LAMMPS dump.Final")
    print(f"{'='*60}\n")
    
    # Parse input
    config = parse_cutpos_input()
    
    if config['n_layers'] is None:
        print("[FAIL] n_layers not specified in cutpos.inp!")
        sys.exit(1)
    
    n_layers = config['n_layers']
    dump_file = config['lammps_dump']
    orthocell = config['orthocell_12atom_sw']
    
    print(f"Configuration:")
    print(f"  n_layers: {n_layers}")
    print(f"  lammps_dump: {dump_file}")
    print(f"  orthocell_12atom_sw: {orthocell}")
    
    # Read dump file
    print(f"\nReading LAMMPS dump: {dump_file}")
    try:
        atoms, lammps_types = parse_lammps_dump(dump_file)
        print(f"  Read {len(atoms)} atoms")
        print(f"  Extracted {len(lammps_types)} LAMMPS types")
        print(f"  Unique types: {sorted(set(lammps_types))}")
    except Exception as e:
        print(f"[FAIL] Error reading dump file: {e}")
        sys.exit(1)
    
    # Load layer definitions
    print(f"\nLoading layer definitions...")
    layers = []
    for ilayer in range(n_layers):
        layer_file = f'layer{ilayer+1}.inp'
        try:
            layer = Layer(layer_file)
            layers.append(layer)
            print(f"  Loaded {layer_file}")
        except Exception as e:
            print(f"[FAIL] Error loading {layer_file}: {e}")
            sys.exit(1)
    
    # Build tag-to-layer mapping
    print(f"\nBuilding tag-to-layer mapping...")
    type_to_symbol = {}
    type_to_layer = {}
    layer_tags = {i: set() for i in range(n_layers)}
    
    for ilayer, layer in enumerate(layers):
        if orthocell and layer.has_orthocell:
            uc_symbols = layer.unitcell_ortho.get_chemical_symbols()
            uc_tags = layer.unitcell_ortho.get_tags()
        else:
            uc_symbols = layer.unitcell.get_chemical_symbols()
            uc_tags = layer.unitcell.get_tags()
        
        for tag, sym in zip(uc_tags, uc_symbols):
            tag_int = int(tag)
            if tag_int in type_to_layer:
                print(f"[FAIL] Tag {tag_int} appears in multiple layers! Make sure that the basis atoms in the layer input files have unique tags.")
                sys.exit(1)
            type_to_symbol[tag_int] = sym
            type_to_layer[tag_int] = ilayer
            layer_tags[ilayer].add(tag_int)
    
    for ilayer in range(n_layers):
        tags_list = sorted(layer_tags[ilayer])
        symbols_in_layer = set([type_to_symbol[t] for t in tags_list])
        print(f"  Layer {ilayer+1}: tags {tags_list} -> {sorted(symbols_in_layer)}")
    
    # Validate tags
    expected_tags = set(type_to_symbol.keys())
    lammps_tags_set = set(lammps_types)
    if expected_tags != lammps_tags_set:
        print(f"[FAIL] Tag mismatch!")
        print(f"  Expected: {sorted(expected_tags)}")
        print(f"  Found: {sorted(lammps_tags_set)}")
        sys.exit(1)
    print(f"  Tag validation passed")
    
    # Set chemical symbols
    symbols = [type_to_symbol[atype] for atype in lammps_types]
    atoms.set_chemical_symbols(symbols)
    atoms.set_tags(lammps_types)
    
    # Write full relaxed structure
    write('relaxed_structure.cif', atoms)
    print(f"\nWrote relaxed_structure.cif ({len(atoms)} atoms)")
    
    # Cut to smaller supercell if specified
    if config['cut_vectors'] is not None and config['lattice_parameters'] is not None:
        print(f"\nCutting to smaller supercell...")
        cut_vectors_cart = config['cut_vectors'] * config['lattice_parameters']
        struct_cut = hexcut(atoms, cut_vectors_cart, atoms.get_cell())
        write('cut_structure.cif', struct_cut)
        print(f"  Wrote cut_structure.cif ({len(struct_cut)} atoms)")
    else:
        struct_cut = atoms
        print(f"\nNo cut specified, using full structure")
    
    # Extract individual layers
    print(f"\nExtracting individual layers...")
    cut_positions = struct_cut.get_positions()
    cut_symbols = struct_cut.get_chemical_symbols()
    cut_types = struct_cut.get_tags()
    cut_cell = struct_cut.get_cell()
    
    for ilayer in range(n_layers):
        # Find atoms belonging to this layer
        layer_mask = np.array([atype in layer_tags[ilayer] for atype in cut_types])
        layer_indices = np.where(layer_mask)[0]
        
        if len(layer_indices) == 0:
            print(f"  WARNING: Layer {ilayer+1}: No atoms found")
            continue
        
        # Create ASE Atoms object for this layer
        layer_positions = cut_positions[layer_indices]
        layer_symbols_list = [cut_symbols[i] for i in layer_indices]
        
        layer_atoms = Atoms(
            symbols=layer_symbols_list,
            positions=layer_positions,
            cell=cut_cell,
            pbc=[True, True, True]
        )
        
        # Write layer CIF file
        output_cif = f'layer_{ilayer+1}.cif'
        write(output_cif, layer_atoms)
        print(f"  Layer {ilayer+1}: {len(layer_indices)} atoms -> {output_cif}")
        
        # Also write coordinate file
        output_coords = f'layer_{ilayer+1}_coords.dat'
        with open(output_coords, 'w') as f:
            for idx in layer_indices:
                pos = cut_positions[idx]
                sym = cut_symbols[idx]
                f.write(f"{sym} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n")
        print(f"             {len(layer_indices)} atoms -> {output_coords}")
    
    print(f"\n{'='*60}")
    print("COMPLETE: Extracted all layers")
    print(f"{'='*60}")
    print(f"\nNext step: Run analysis workflow")
    print(f"  python ../../Src/run_analysis.py")


if __name__ == '__main__':
    main()
