import numpy as np
from .material_detector import detect_material_type


def generate_lammps_input(structure_file, layers, config, use_orthocell=False):
    """
    Generate lammps.in file based on materials present in the system.
    
    Args:
        structure_file: Path to LAMMPS structure data file
        layers: List of Layer objects
        config: Configuration dict from twisterase.inp
        use_orthocell: Whether orthorhombic cells are used
    """
    # Detect materials and gather information
    material_info = _analyze_system(structure_file, layers, use_orthocell)
    
    # Determine which generator to use based on materials present
    has_hbn = material_info['has_hbn']
    has_graphene = material_info['has_graphene']
    has_tmd = material_info['has_tmd']
    
    if has_hbn and not has_graphene and not has_tmd:
        # Pure hBN system
        _generate_hbn_input(structure_file, material_info, config)
    elif has_graphene and not has_hbn and not has_tmd:
        # Pure graphene system
        _generate_graphene_input(structure_file, material_info, config)
    elif has_tmd and not has_hbn and not has_graphene:
        # Pure TMD system
        _generate_tmd_input(structure_file, material_info, config, use_orthocell)
    elif has_hbn or (has_graphene and has_tmd):
        # Mixed system with hBN, or graphene+TMD
        _generate_mixed_input(structure_file, material_info, config)
    else:
        raise ValueError(f"Unsupported material combination: {material_info['layer_types']}")
    
    print(f"Written lammps.in for materials: {', '.join(material_info['layer_types'])}")


def _get_supercell(layer, use_orthocell):
    """Return the appropriate supercell (ortho or hex) for a layer."""
    return layer.supercell_ortho if (use_orthocell and layer.has_orthocell) else layer.supercell


def _get_input_tags(layer, use_orthocell):
    """Return input tags from layer data, preferring ortho tags when appropriate."""
    if use_orthocell and layer.has_orthocell and layer.layer_data.get('atom_tags_ortho'):
        return layer.layer_data.get('atom_tags_ortho')
    return layer.layer_data.get('atom_tags')


def _analyze_system(structure_file, layers, use_orthocell):
    """Analyze the system and gather all necessary information."""
    layer_info = []
    all_symbols = []
    
    # First pass: find max TMD type from layer input files
    max_tmd_type = 0
    for layer in layers:
        supercell = _get_supercell(layer, use_orthocell)
        material_type = detect_material_type(supercell.get_chemical_symbols())
        
        if material_type != 'hBN':
            input_tags = _get_input_tags(layer, use_orthocell)
            if input_tags:
                max_tmd_type = max(max_tmd_type, max(input_tags))
    
    # Second pass: assign atom types
    hbn_layer_count = 0
    for i, layer in enumerate(layers):
        supercell = _get_supercell(layer, use_orthocell)
        symbols = supercell.get_chemical_symbols()
        tags = supercell.get_tags()
        material_type = detect_material_type(symbols)
        
        if material_type == 'hBN':
            input_tags = _get_input_tags(layer, use_orthocell)
            if input_tags:
                unique_tags = sorted(set(input_tags))
            else:
                type_B = max_tmd_type + 2 * hbn_layer_count + 1
                type_N = max_tmd_type + 2 * hbn_layer_count + 2
                unique_tags = [type_B, type_N]
            hbn_layer_count += 1
        else:
            input_tags = _get_input_tags(layer, use_orthocell)
            unique_tags = list(set(input_tags)) if input_tags else list(set(tags))
        
        layer_info.append({
            'index': i,
            'material_type': material_type,
            'symbols': symbols,
            'unique_symbols': list(set(symbols)),
            'tags': tags,
            'unique_tags': unique_tags,
            'num_atoms': len(symbols),
            'intralayer_potential': layer.get_intralayer_potential()
        })
        all_symbols.extend(symbols)
    
    # Read structure file to get atom types and masses
    atom_types, masses, type_to_symbol = _read_structure_info(structure_file)
    layer_types = [info['material_type'] for info in layer_info]
    
    return {
        'n_layers': len(layers),
        'layer_info': layer_info,
        'layer_types': layer_types,
        'has_hbn': 'hBN' in layer_types,
        'has_graphene': 'Graphene' in layer_types,
        'has_tmd': 'TMD' in layer_types,
        'all_symbols': all_symbols,
        'unique_symbols': list(set(all_symbols)),
        'atom_types': atom_types,
        'masses': masses,
        'type_to_symbol': type_to_symbol
    }


def _read_structure_info(structure_file):
    """Read atom types and masses from structure file."""
    atom_types = []
    masses = []
    type_to_symbol = {}
    
    with open(structure_file, 'r') as f:
        lines = f.readlines()
    
    read_masses = False
    for line in lines:
        line = line.strip()
        if "Masses" in line:
            read_masses = True
            continue
        if "Atoms" in line:
            break
        if read_masses and line:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    atype = int(parts[0])
                    mass = float(parts[1])
                    atom_types.append(atype)
                    masses.append(mass)
                    
                    # Map mass to symbol
                    mass_map = {
                        183.84: 'W', 95.95: 'Mo', 78.971: 'Se', 32.06: 'S',
                        127.60: 'Te', 12.011: 'C', 10.81: 'B', 14.007: 'N'
                    }
                    closest = min(mass_map.keys(), key=lambda x: abs(x - mass))
                    type_to_symbol[atype] = mass_map[closest]
                except ValueError:
                    continue
    
    return atom_types, masses, type_to_symbol


def _generate_hbn_input(structure_file, material_info, config):
    """Generate lammps.in for pure hBN system."""
    with open('lammps.in', 'w') as f:
        f.write("# Initialization\n")
        f.write("units           metal\n")
        f.write("dimension       3\n")
        f.write("box tilt        large\n")
        f.write("boundary        p p p\n")
        f.write("atom_style      full\n\n")
        
        f.write("# System and atom definition\n")
        f.write(f"read_data       {structure_file}\n\n")
        
        for atype, mass in zip(material_info['atom_types'], material_info['masses']):
            f.write(f"mass            {atype} {mass}\n")
        f.write("\n")
        
        # Group definitions
        for i in range(material_info['n_layers']):
            f.write(f"group           layer{i+1} molecule {i+1}\n")
        f.write("\n")
        
        f.write("# Potentials\n")
        f.write("pair_style      hybrid/overlay tersoff ilp/graphene/hbn 16.0 coul/shield 16.0 1\n\n")
        
        f.write("# Intralayer Interaction\n")
        pot = material_info['layer_info'][0]['intralayer_potential']
        symbols_for_types = [material_info['type_to_symbol'][atype] for atype in material_info['atom_types']]
        f.write(f"pair_coeff      * * tersoff {pot} {' '.join(symbols_for_types)}\n\n")
        
        f.write("# Interlayer Interaction\n")
        f.write(f"pair_coeff      * * ilp/graphene/hbn BNCH.ILP {' '.join(symbols_for_types)}\n\n")
        
        # Coulomb interactions
        f.write("# Coulomb interactions between layers\n")
        for i in range(material_info['n_layers'] - 1):
            for j in range(i+1, material_info['n_layers']):
                layer_i_types = material_info['layer_info'][i]['unique_tags']
                layer_j_types = material_info['layer_info'][j]['unique_tags']
                for type_i in layer_i_types:
                    for type_j in layer_j_types:
                        sym_i = material_info['type_to_symbol'][type_i]
                        sym_j = material_info['type_to_symbol'][type_j]
                        if sym_i == 'B' and sym_j == 'B':
                            coeff = 0.70
                        elif (sym_i == 'B' and sym_j == 'N') or (sym_i == 'N' and sym_j == 'B'):
                            coeff = 0.69498201415576216335
                        elif sym_i == 'N' and sym_j == 'N':
                            coeff = 0.69
                        else:
                            coeff = 0.0
                        f.write(f"pair_coeff      {type_i} {type_j} coul/shield {coeff}\n")
        f.write("\n")
        
        f.write("# Neighbor settings\n")
        f.write("neighbor        2.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")
        
        f.write("# Optimization\n")
        f.write("dump            1 all custom 100 dump.minimization id type x y z\n")
        f.write("thermo          1000\n")
        f.write("thermo_style    custom step pe press\n")
        if material_info["has_graphene"]:
            f.write("min_style       cg\n")
            f.write("minimize        1.0e-10 1.0e-12 100000 1000000\n")
        else:
            f.write("min_style       fire\n")
            f.write("minimize        1.0e-10 1.0e-12 100000 1000000\n")
        f.write("undump          1\n")
        f.write("write_dump      all custom dump.Final id type x y z modify sort id\n")
        f.write('print "Done!"\n')


def _generate_graphene_input(structure_file, material_info, config):
    """Generate lammps.in for pure graphene system."""
    with open('lammps.in', 'w') as f:
        f.write("#1 General\n")
        f.write("units           metal\n")
        f.write("dimension       3\n")
        f.write("atom_style      atomic\n")
        f.write("neighbor        0.3 bin\n\n")
        
        f.write("#2 Structure\n")
        f.write("boundary        p p p\n")
        f.write(f"read_data       {structure_file}\n")
        for atype, mass in zip(material_info['atom_types'], material_info['masses']):
            f.write(f"mass {atype} {mass}\n")
        f.write("\n")
        
        f.write("#3 Potential definitions\n")
        f.write("pair_style      hybrid/overlay rebo kolmogorov/crespi/z 14.0\n\n")
        
        f.write("#4 Intralayer interactions\n")
        pot = material_info['layer_info'][0]['intralayer_potential']
        # REBO requires one element per atom type
        symbols_for_types = [material_info['type_to_symbol'][atype] for atype in material_info['atom_types']]
        f.write(f"pair_coeff      * * rebo {pot} {' '.join(symbols_for_types)}\n\n")
        
        f.write("#5 Interlayer interactions\n")
        inter_pot = config.get('interlayer_potential', 'CC.KC')
        symbols_for_types = [material_info['type_to_symbol'][atype] for atype in material_info['atom_types']]
        for i in range(material_info['n_layers'] - 1):
            type_i = sorted(material_info['layer_info'][i]['unique_tags'])[0]
            type_j = sorted(material_info['layer_info'][i+1]['unique_tags'])[0]
            f.write(f"pair_coeff      {type_i} {type_j} kolmogorov/crespi/z {inter_pot} {' '.join(symbols_for_types)}\n")
        f.write("\n")
        
        f.write("#6 Optimize at 0 K\n")
        f.write("dump            1 all custom 400 dump.minimization id type x y z\n")
        if material_info["has_graphene"]:
            f.write("min_style       cg\n")
            f.write("minimize        1.0e-10 1.0e-12 100000 1000000\n")
        else:
            f.write("min_style       fire\n")
            f.write("minimize        0.0 1.0e-8 1000000 1000000\n")
        f.write("undump          1\n")
        f.write("write_dump      all custom dump.Final id type x y z modify sort id\n")
        f.write('print "Done!"\n')


def _generate_tmd_input(structure_file, material_info, config, use_orthocell):
    """Generate lammps.in for pure TMD system (handles both hex and ortho basis)."""
    n_layers = material_info['n_layers']
    
    # Determine basis type from number of unique atom types per layer
    # Hexagonal: typically 3 types (1 metal + 2 chalcogens)
    # Orthorhombic: typically 12 types (4 metals + 8 chalcogens)
    basis_size = len(material_info['layer_info'][0]['unique_tags'])
    is_ortho_basis = basis_size >= 12
    
    if is_ortho_basis:
        print(f"Detected orthorhombic TMD basis ({basis_size} atom types per layer)")
        # For ortho: 64 interactions per layer pair (4 sets × 4×4)
        num_interactions = 64 * (n_layers - 1)
    else:
        print(f"Detected hexagonal TMD basis ({basis_size} atom types per layer)")
        # For hex: classify and count interactions
        layer_classifications = []
        for i, layer_info in enumerate(material_info['layer_info']):
            classifications = _classify_tmd_atoms(layer_info, structure_file)
            layer_classifications.append(classifications)
        num_interactions = _count_tmd_interactions(layer_classifications, material_info)
    
    print(f"  Total interlayer interactions: {num_interactions}")
    
    with open('lammps.in', 'w') as f:
        f.write("#1 General\n")
        f.write("units           metal\n")
        f.write("dimension       3\n")
        f.write("atom_style      atomic\n")
        f.write("neighbor        0.3 bin\n\n")
        
        f.write("#2 Structure\n")
        f.write("boundary        p p p\n")
        f.write(f"read_data       {structure_file}\n")
        for atype, mass in zip(material_info['atom_types'], material_info['masses']):
            f.write(f"mass {atype} {mass}\n")
        f.write("\n")
        
        f.write("#3 Potential definitions\n")
        pair_style = f"pair_style hybrid/overlay {' '.join(['sw/mod'] * n_layers)} "
        pair_style += ' '.join(['kolmogorov/crespi/z 14.0'] * num_interactions) + " lj/cut 10.0\n"
        f.write(pair_style)
        f.write("\n")
        
        f.write("#4 Intralayer interactions\n")
        for i, layer_info in enumerate(material_info['layer_info']):
            symbols_line = ['NULL'] * len(material_info['atom_types'])
            
            if is_ortho_basis:
                # Orthorhombic: use numbered symbols (W1, Se1, Se2, etc.)
                symbol_counts = {}
                sorted_tags = sorted(layer_info['unique_tags'])
                
                for tag in sorted_tags:
                    base_symbol = material_info['type_to_symbol'][tag]
                    if base_symbol not in symbol_counts:
                        symbol_counts[base_symbol] = 0
                    symbol_counts[base_symbol] += 1
                    numbered_symbol = f"{base_symbol}{symbol_counts[base_symbol]}"
                    symbols_line[tag - 1] = numbered_symbol
            else:
                # Hexagonal: use plain symbols (W, Se, Se)
                for tag in layer_info['unique_tags']:
                    symbols_line[tag - 1] = material_info['type_to_symbol'][tag]
            
            pot = layer_info['intralayer_potential']
            f.write(f"pair_coeff * * sw/mod {i+1} {pot} {' '.join(symbols_line)}\n")
        f.write("\n")
        
        f.write("#5 Interlayer interactions\n")
        inter_pot = config.get('interlayer_potential', 'TMD.KC')
        interaction_num = 1
        
        if is_ortho_basis:
            # Orthorhombic: Generate 64 interactions per layer pair
            _generate_ortho_interlayer_interactions(f, material_info, inter_pot, interaction_num, structure_file)
        else:
            # Hexagonal: Use classification-based approach
            for i in range(n_layers - 1):
                for j, (tag_i, class_i, sym_i) in enumerate(layer_classifications[i]):
                    for k, (tag_j, class_j, sym_j) in enumerate(layer_classifications[i+1]):
                        if _should_interact_tmd(class_i, class_j):
                            symbols_line = ['NULL'] * len(material_info['atom_types'])
                            symbols_line[tag_i - 1] = sym_i
                            symbols_line[tag_j - 1] = sym_j
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
        f.write("pair_coeff * * lj/cut 0.0 3.4\n\n")
        
        f.write("#6 Optimize at 0 K\n")
        f.write("dump            1 all custom 400 dump.minimization id type x y z\n")
        f.write("min_style       fire\n")
        f.write("minimize        0.0 1.0e-8 1000000 1000000\n")
        f.write("undump          1\n")
        f.write("write_dump      all custom dump.Final id type x y z modify sort id\n")
        f.write('print "Done!"\n')


def _generate_mixed_input(structure_file, material_info, config):
    """Generate lammps.in for mixed material systems (hBN + TMD/Graphene)."""
    if material_info['has_hbn'] and material_info['has_tmd']:
        from .input_generator_hbn_tmd import generate_hbn_tmd_input
        generate_hbn_tmd_input(structure_file, material_info, config)
    elif material_info['has_hbn'] and material_info['has_graphene']:
        from .input_generator_hbn_tmd import generate_hbn_graphene_input
        generate_hbn_graphene_input(structure_file, material_info, config)
    elif material_info['has_hbn']:
        # Pure hBN system routed through mixed
        _generate_hbn_input(structure_file, material_info, config)
    else:
        raise NotImplementedError("Mixed TMD+Graphene systems not yet implemented")


def _partition_chalcogens_ortho(layer_info, type_to_symbol, structure_file=None):
    """
    Partition chalcogens in orthorhombic structure into two groups based on z-coordinates.
    Group A: chalcogens at lower z-coordinate (X_l)
    Group B: chalcogens at upper z-coordinate (X_u)
    """
    TMD_CHALCOGENS = {'S', 'Se', 'Te'}
    
    # If structure_file provided, use z-coordinate based classification
    if structure_file:
        # Use the same z-based classification as hexagonal TMD
        classifications = _classify_tmd_atoms(layer_info, structure_file)
        
        group_A = []  # X_l (lower chalcogens)
        group_B = []  # X_u (upper chalcogens)
        
        for tag, class_i, sym in classifications:
            if class_i == 'X_l':
                group_A.append(tag)
            elif class_i == 'X_u':
                group_B.append(tag)
        
        return group_A, group_B
    else:
        # Fallback to tag-based approach if no structure file
        chalc_tags = []
        for tag in layer_info['unique_tags']:
            if tag in type_to_symbol and type_to_symbol[tag] in TMD_CHALCOGENS:
                chalc_tags.append(tag)
        
        chalc_tags_sorted = sorted(chalc_tags)
        group_A = chalc_tags_sorted[0::2]
        group_B = chalc_tags_sorted[1::2]
        
        return group_A, group_B


def _partition_metals_ortho(layer_info, type_to_symbol):
    """Get all metal atom tags in orthorhombic structure."""
    TMD_METALS = {'Mo', 'W', 'Nb', 'Ta', 'Re', 'Tc'}
    
    # Find metal tags using type_to_symbol mapping
    metal_tags = []
    for tag in layer_info['unique_tags']:
        if tag in type_to_symbol and type_to_symbol[tag] in TMD_METALS:
            metal_tags.append(tag)
    
    return sorted(metal_tags)


def _write_kc_pairs(f, tags_i, tags_j, type_to_symbol, n_types, inter_pot, interaction_num, sym_override=None):
    """
    Write KC pair_coeff lines for all (tag_i, tag_j) pairs.
    
    Args:
        sym_override: dict mapping tag -> symbol override (e.g. for hBN atoms mapped to 'B')
    Returns:
        Updated interaction_num.
    """
    if sym_override is None:
        sym_override = {}
    for tag_i in tags_i:
        for tag_j in tags_j:
            symbols_line = ['NULL'] * n_types
            symbols_line[tag_i - 1] = sym_override.get(tag_i, type_to_symbol[tag_i])
            symbols_line[tag_j - 1] = sym_override.get(tag_j, type_to_symbol[tag_j])
            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
            interaction_num += 1
    return interaction_num


def _generate_ortho_interlayer_interactions(f, material_info, inter_pot, start_interaction_num, structure_file):
    """
    Generate 64 interlayer interactions per layer pair for orthorhombic TMD structures.
    
    Following the original TwisterASE approach with z-based chalcogen partitioning:
    - Set 1: Lower chalc_B (X_u) × Upper chalc_A (X_l) (4×4 = 16)
    - Set 2: Lower metals × Upper metals (4×4 = 16)
    - Set 3: Lower chalc_B (X_u) × Upper metals (4×4 = 16)
    - Set 4: Lower metals × Upper chalc_A (X_l) (4×4 = 16)
    Total: 64 interactions per layer pair
    """
    n_layers = material_info['n_layers']
    type_to_symbol = material_info['type_to_symbol']
    n_types = len(material_info['atom_types'])
    interaction_num = start_interaction_num
    
    for i in range(n_layers - 1):
        lower_layer = material_info['layer_info'][i]
        upper_layer = material_info['layer_info'][i + 1]
        
        lower_chalc_A, lower_chalc_B = _partition_chalcogens_ortho(lower_layer, type_to_symbol, structure_file)
        upper_chalc_A, upper_chalc_B = _partition_chalcogens_ortho(upper_layer, type_to_symbol, structure_file)
        lower_metals = _partition_metals_ortho(lower_layer, type_to_symbol)
        upper_metals = _partition_metals_ortho(upper_layer, type_to_symbol)
        
        interaction_num = _write_kc_pairs(f, lower_chalc_B, upper_chalc_A, type_to_symbol, n_types, inter_pot, interaction_num)
        interaction_num = _write_kc_pairs(f, lower_metals, upper_metals, type_to_symbol, n_types, inter_pot, interaction_num)
        interaction_num = _write_kc_pairs(f, lower_chalc_B, upper_metals, type_to_symbol, n_types, inter_pot, interaction_num)
        interaction_num = _write_kc_pairs(f, lower_metals, upper_chalc_A, type_to_symbol, n_types, inter_pot, interaction_num)


def _classify_tmd_atoms(layer_info, structure_file):
    """
    Classify TMD atoms as TM, X_u, or X_l based on z-coordinates.
    
    X_u: Chalcogen with maximum z-coordinate (upper)
    X_l: Chalcogen with lower z-coordinate (lower)
    TM: Transition metal
    """
    TMD_METALS = {'Mo', 'W', 'Nb', 'Ta', 'Re', 'Tc'}
    TMD_CHALCOGENS = {'S', 'Se', 'Te'}
    
    # Read atom positions from structure file to get z-coordinates
    # Store ALL occurrences of each atom type
    atom_z_coords = {}
    with open(structure_file, 'r') as f:
        lines = f.readlines()
    
    read_atoms = False
    for line in lines:
        line = line.strip()
        if "Atoms" in line:
            read_atoms = True
            continue
        if read_atoms and line:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    atom_id = int(parts[0])
                    atom_type = int(parts[1])
                    # Handle both atomic and full formats
                    if len(parts) == 5:  # atomic format: id type x y z
                        z = float(parts[4])
                    else:  # full format: id mol type charge x y z
                        z = float(parts[6])
                    
                    # Store all z-coordinates for this type
                    if atom_type not in atom_z_coords:
                        atom_z_coords[atom_type] = []
                    atom_z_coords[atom_type].append(z)
                except (ValueError, IndexError):
                    continue
    
    # Average z-coordinates for each atom type
    atom_z_avg = {}
    for atype, z_list in atom_z_coords.items():
        atom_z_avg[atype] = np.mean(z_list)
    
    # Read type_to_symbol from structure file
    type_to_symbol = {}
    atom_types, masses, type_to_symbol = _read_structure_info(structure_file)
    
    # Classify atoms based on symbol and z-coordinate
    classifications = []
    
    # Group chalcogens by z-coordinate to find max
    chalcogen_tags = []
    chalcogen_z = []
    for tag in layer_info['unique_tags']:
        if tag in type_to_symbol:
            sym = type_to_symbol[tag]
            if sym in TMD_CHALCOGENS and tag in atom_z_avg:
                chalcogen_tags.append(tag)
                chalcogen_z.append(atom_z_avg[tag])
    
    # Find maximum z among chalcogens
    max_z = max(chalcogen_z) if chalcogen_z else 0
    
    # Classify all atoms
    for tag in layer_info['unique_tags']:
        if tag in type_to_symbol:
            sym = type_to_symbol[tag]
            if sym in TMD_METALS:
                classifications.append((tag, 'TM', sym))
            elif sym in TMD_CHALCOGENS:
                if tag in atom_z_avg:
                    z = atom_z_avg[tag]
                    # X_u if at maximum z, X_l otherwise
                    if abs(z - max_z) < 0.01:  # Increased tolerance for ortho structures
                        classifications.append((tag, 'X_u', sym))
                    else:
                        classifications.append((tag, 'X_l', sym))
                else:
                    classifications.append((tag, 'X', sym))
            else:
                classifications.append((tag, 'X', sym))
    
    return classifications


def _should_interact_tmd(class1, class2):
    """Determine if two TMD atom types should have interlayer interaction."""
    return ((class1 == 'TM' and class2 == 'TM') or
            (class1 == 'TM' and class2 == 'X_l') or
            (class1 == 'X_u' and class2 in ['TM', 'X_l']))


def _count_tmd_interactions(layer_classifications, material_info):
    """Count number of TMD interlayer interactions."""
    count = 0
    for i in range(len(layer_classifications) - 1):
        for tag_i, class_i, sym_i in layer_classifications[i]:
            for tag_j, class_j, sym_j in layer_classifications[i+1]:
                if _should_interact_tmd(class_i, class_j):
                    count += 1
    return count
