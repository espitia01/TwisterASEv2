import numpy as np
from ase import Atoms
from ase.io import write
from .material_detector import detect_system_materials, detect_material_type


def restricted_triclinic(lattice_vectors):
    """Convert lattice vectors to LAMMPS restricted triclinic format."""
    a, b, c = np.array(lattice_vectors[0]), np.array(lattice_vectors[1]), np.array(lattice_vectors[2])
    
    e1 = a / np.linalg.norm(a)
    
    b_proj = np.dot(b, e1) * e1
    b_perp = b - b_proj
    if np.linalg.norm(b_perp) < 1e-10:
        raise ValueError("a and b are colinear; no unique basis")
    e2 = b_perp / np.linalg.norm(b_perp)
    
    e3 = np.cross(e1, e2)
    
    if np.dot(e3, c) <= 0:
        raise ValueError("c' would have non-positive z-component")
    
    O = np.vstack([e1, e2, e3])
    return (O @ lattice_vectors.T).T


def write_structure_file(filename, layers, use_orthocell=False):
    """
    Write LAMMPS structure file with appropriate format based on materials present.
    
    Args:
        filename: Output filename
        layers: List of Layer objects
        use_orthocell: If True, use orthorhombic supercells
    """
    # Detect material types
    material_info = detect_system_materials(layers)
    
    # Combine layers into single structure
    if use_orthocell:
        first = True
        combined = None
        for layer in layers:
            sc = layer.supercell_ortho if layer.supercell_ortho is not None else layer.supercell
            if first:
                combined = sc.copy()
                first = False
            else:
                combined += sc
    else:
        combined = layers[0].supercell.copy()
        for i in range(1, len(layers)):
            combined += layers[i].supercell
    
    # Convert to restricted triclinic
    old_lattice_vectors = np.array(combined.get_cell())
    lattice_vectors = restricted_triclinic(old_lattice_vectors)
    
    # Determine format: use atomic for mixed hBN+TMD with ortho hBN, full for pure hBN
    has_hbn = material_info['has_hbn']
    has_tmd = material_info['has_tmd'] or material_info['has_graphene']
    is_mixed_ortho_hbn = has_hbn and has_tmd and use_orthocell
    
    # Write appropriate format
    if has_hbn and not has_tmd:
        # Pure hBN: use full format with molecule ID and charge
        pos_new = _write_full_format(filename, layers, combined, lattice_vectors, use_orthocell)
    elif is_mixed_ortho_hbn:
        # Mixed hBN+TMD with ortho: use atomic format
        # Set correct tags on combined object from layer ortho tags
        combined = _set_combined_tags(combined, layers, use_orthocell)
        pos_new = _write_atomic_format(filename, combined, lattice_vectors)
    elif material_info['requires_full_format']:
        # hBN present but not ortho: use full format
        pos_new = _write_full_format(filename, layers, combined, lattice_vectors, use_orthocell)
    else:
        # TMD/Graphene only: use atomic format
        pos_new = _write_atomic_format(filename, combined, lattice_vectors)
    
    # Write CIF file for visualization
    struct_new = Atoms(
        symbols=combined.get_chemical_symbols(),
        positions=pos_new,
        cell=lattice_vectors,
        pbc=[True, True, True]
    )
    write('superlattice_lammps.cif', struct_new)
    print(f"Written superlattice_lammps.cif with {len(struct_new)} atoms")
    
    return material_info


def _set_combined_tags(combined, layers, use_orthocell):
    """
    Set correct atom tags on the combined Atoms object for mixed hBN+TMD systems.
    make_supercell preserves tags from the unitcell for all material types.
    """
    all_tags = []
    for layer in layers:
        layer_atoms = layer.supercell_ortho if (use_orthocell and layer.has_orthocell) else layer.supercell
        all_tags.extend(layer_atoms.get_tags().tolist())
    
    combined.set_tags(all_tags)
    return combined


def _write_atomic_format(filename, atoms, lattice_vectors):
    """
    Write LAMMPS data file in atomic format (for TMD/Graphene).
    Format: atom-ID atom-type x y z
    """
    with open(filename, 'w') as f:
        f.write('LAMMPS data file\n\n')
        f.write('{} atoms\n'.format(len(atoms.positions)))
        f.write('{} atom types\n'.format(len(np.unique(atoms.get_tags()))))
        
        # Box bounds
        f.write('0.0 {} xlo xhi\n'.format(lattice_vectors[0, 0]))
        f.write('0.0 {} ylo yhi\n'.format(lattice_vectors[1, 1]))
        f.write('0.0 {} zlo zhi\n'.format(lattice_vectors[2, 2]))
        f.write('{} {} {} xy xz yz\n'.format(lattice_vectors[1, 0], lattice_vectors[2, 0], lattice_vectors[2, 1]))
        f.write('\n')
        
        # Masses
        f.write('Masses\n\n')
        unique_tags, indices = np.unique(atoms.get_tags(), return_index=True)
        for ind in indices:
            f.write('{} {}\n'.format(atoms.get_tags()[ind], atoms.get_masses()[ind]))
        f.write('\n')
        
        # Atoms
        f.write('Atoms\n\n')
        tags = atoms.get_tags()
        pos_arr = np.zeros((len(atoms.positions), 3))
        for i, iat in enumerate(atoms.get_scaled_positions()):
            pos = np.dot(iat, lattice_vectors)
            pos_arr[i] = pos
            f.write('{} {} {} {} {}\n'.format(i + 1, tags[i], pos[0], pos[1], pos[2]))
    
    return pos_arr


def _write_full_format(filename, layers, combined, lattice_vectors, use_orthocell):
    """
    Write LAMMPS data file in full format (for systems with hBN).
    Format: atom-ID molecule-ID atom-type charge x y z
    """
    # Assign molecule IDs (layer numbers)
    molecule_ids = []
    for i, layer in enumerate(layers):
        if use_orthocell and layer.has_orthocell:
            n_atoms = len(layer.supercell_ortho)
        else:
            n_atoms = len(layer.supercell)
        molecule_ids.extend([i + 1] * n_atoms)
    molecule_ids = np.array(molecule_ids)
    
    # Get chemical symbols
    chem_syms = combined.get_chemical_symbols()
    
    # Compute atom types and charges based on material type and layer
    # For mixed systems, assign TMD types first, then hBN types after
    atom_types = []
    charges = []
    
    def _get_sc(layer):
        return layer.supercell_ortho if (use_orthocell and layer.has_orthocell) else layer.supercell
    
    def _get_tags(layer):
        if use_orthocell and layer.has_orthocell and layer.layer_data.get('atom_tags_ortho'):
            return layer.layer_data.get('atom_tags_ortho')
        return layer.layer_data.get('atom_tags')
    
    # First pass: determine max TMD type
    max_tmd_type = 0
    for layer in layers:
        layer_atoms = _get_sc(layer)
        material_type = detect_material_type(layer_atoms.get_chemical_symbols())
        
        if material_type != 'hBN':
            input_tags = _get_tags(layer)
            if input_tags:
                max_tmd_type = max(max_tmd_type, max(input_tags))
    
    # Second pass: assign atom types
    hbn_layer_count = 0
    for layer_idx, layer in enumerate(layers):
        layer_atoms = _get_sc(layer)
        layer_symbols = layer_atoms.get_chemical_symbols()
        material_type = detect_material_type(layer_symbols)
        
        if material_type == 'hBN':
            # Check if hBN has explicit ortho tags
            if use_orthocell and layer.has_orthocell and layer.layer_data.get('atom_tags_ortho'):
                # Use ortho tags from input file - supercell preserves them via make_supercell
                layer_tags = layer_atoms.get_tags()
                # Build symbol-to-charge map
                for i_atom, sym in enumerate(layer_symbols):
                    atom_types.append(int(layer_tags[i_atom]))
                    if sym == 'B':
                        charges.append(0.42)
                    elif sym == 'N':
                        charges.append(-0.42)
                    else:
                        charges.append(0.0)
                hbn_layer_count += 1
            else:
                # Fallback: assign types after all TMD types
                # B gets odd type, N gets even type
                type_B = max_tmd_type + 2 * hbn_layer_count + 1
                type_N = max_tmd_type + 2 * hbn_layer_count + 2
                hbn_layer_count += 1
                
                for sym in layer_symbols:
                    if sym == 'B':
                        atom_types.append(type_B)
                        charges.append(0.42)
                    elif sym == 'N':
                        atom_types.append(type_N)
                        charges.append(-0.42)
        else:
            # TMD/Graphene: use original tags from layer input file
            if use_orthocell and layer.has_orthocell and layer.layer_data.get('atom_tags_ortho'):
                input_tags = layer.layer_data.get('atom_tags_ortho')
                input_symbols = layer.layer_data.get('atom_symbols_ortho')
                input_positions = layer.layer_data.get('atom_positions_ortho')
            else:
                input_tags = layer.layer_data.get('atom_tags')
                input_symbols = layer.layer_data.get('atom_symbols')
                input_positions = layer.layer_data.get('atom_positions')
            
            if input_tags and input_symbols and input_positions is not None:
                is_ortho_basis = len(set(input_tags)) >= 12
                
                if is_ortho_basis:
                    # Orthorhombic basis: use supercell tags directly
                    # make_supercell preserves tags for ortho cells
                    layer_tags = layer_atoms.get_tags()
                    for tag in layer_tags:
                        atom_types.append(int(tag))
                        charges.append(0.0)
                else:
                    # Hexagonal basis: use z-coordinate classification
                    TMD_METALS = {'Mo', 'W', 'Nb', 'Ta', 'Re', 'Tc'}
                    TMD_CHALCOGENS = {'S', 'Se', 'Te'}
                    
                    # Get z-coordinates from input file (scaled positions)
                    input_z = [pos[2] for pos in input_positions]
                    
                    # Create symbol-to-tag mapping based on z-coordinate
                    metal_tag = None
                    chalc_tags = []  # [(z, tag), ...]
                    
                    for sym, z, tag in zip(input_symbols, input_z, input_tags):
                        if sym in TMD_METALS:
                            metal_tag = tag
                        elif sym in TMD_CHALCOGENS:
                            chalc_tags.append((z, tag))
                    
                    # Sort chalcogens by z to determine lower/upper
                    chalc_tags.sort(key=lambda x: x[0])
                    chalc_lower_tag = chalc_tags[0][1] if len(chalc_tags) > 0 else None
                    chalc_upper_tag = chalc_tags[1][1] if len(chalc_tags) > 1 else chalc_lower_tag
                    
                    # Get supercell positions and compute z-thresholds
                    supercell_positions = layer_atoms.get_positions()
                    supercell_z = supercell_positions[:, 2]
                    
                    metal_z_avg = np.mean([supercell_z[i] for i, sym in enumerate(layer_symbols) if sym in TMD_METALS])
                    chalc_z_list = [supercell_z[i] for i, sym in enumerate(layer_symbols) if sym in TMD_CHALCOGENS]
                    
                    if chalc_z_list:
                        chalc_z_mid = (max(chalc_z_list) + min(chalc_z_list)) / 2
                    else:
                        chalc_z_mid = metal_z_avg
                    
                    # Assign tags based on symbol and z-coordinate
                    for i, sym in enumerate(layer_symbols):
                        z = supercell_z[i]
                        if sym in TMD_METALS:
                            atom_types.append(int(metal_tag))
                        elif sym in TMD_CHALCOGENS:
                            if z < chalc_z_mid:
                                atom_types.append(int(chalc_lower_tag))
                            else:
                                atom_types.append(int(chalc_upper_tag))
                        else:
                            atom_types.append(int(input_tags[0]))
                        charges.append(0.0)
            else:
                # Fallback: use supercell tags
                layer_tags = layer_atoms.get_tags()
                for tag in layer_tags:
                    atom_types.append(int(tag))
                    charges.append(0.0)
    
    atom_types = np.array(atom_types)
    charges = np.array(charges)
    
    # Write file
    with open(filename, 'w') as f:
        f.write('LAMMPS data file\n\n')
        f.write('{} atoms\n'.format(len(combined.positions)))
        f.write('{} atom types\n'.format(len(np.unique(atom_types))))
        
        # Box bounds
        f.write('0.0 {} xlo xhi\n'.format(lattice_vectors[0, 0]))
        f.write('0.0 {} ylo yhi\n'.format(lattice_vectors[1, 1]))
        f.write('0.0 {} zlo zhi\n'.format(lattice_vectors[2, 2]))
        f.write('{} {} {} xy xz yz\n'.format(lattice_vectors[1, 0], lattice_vectors[2, 0], lattice_vectors[2, 1]))
        f.write('\n')
        
        # Masses
        f.write('Masses\n\n')
        unique_types = np.unique(atom_types)
        for atype in sorted(unique_types):
            # Find first atom with this type to get mass
            idx = np.where(atom_types == atype)[0][0]
            mass = combined.get_masses()[idx]
            f.write('{} {}\n'.format(int(atype), mass))
        f.write('\n')
        
        # Atoms (full format)
        f.write('Atoms\n\n')
        pos_arr = np.zeros((len(combined.positions), 3))
        for i, iat in enumerate(combined.get_scaled_positions()):
            pos = np.dot(iat, lattice_vectors)
            pos_arr[i] = pos
            f.write('{} {} {} {} {} {} {}\n'.format(
                i + 1,              # atom ID
                molecule_ids[i],    # molecule ID (layer number)
                int(atom_types[i]), # atom type
                charges[i],         # charge
                pos[0], pos[1], pos[2]
            ))
    
    return pos_arr
