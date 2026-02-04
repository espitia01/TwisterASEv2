"""
Generator for mixed hBN + TMD systems.
Handles heterostructures like hBN-WSe2-WSe2-hBN.
"""

def generate_hbn_tmd_input(structure_file, material_info, config):
    """
    Generate lammps.in for mixed hBN + TMD systems.
    
    Uses:
    - Tersoff for hBN intralayer (separate instance per layer)
    - SW/mod for TMD intralayer (separate instance per layer)
    - KC for all interlayer interactions (hBN-hBN, TMD-TMD, hBN-TMD)
    """
    n_layers = material_info['n_layers']
    
    # Count TMD and hBN layers
    n_hbn = sum(1 for info in material_info['layer_info'] if info['material_type'] == 'hBN')
    n_tmd = sum(1 for info in material_info['layer_info'] if info['material_type'] == 'TMD')
    
    # Count all interlayer interactions
    total_interactions = 0
    TMD_CHALCOGENS = {'S', 'Se', 'Te'}
    
    for i in range(n_layers - 1):
        lower_layer = material_info['layer_info'][i]
        upper_layer = material_info['layer_info'][i+1]
        
        if lower_layer['material_type'] == 'TMD' and upper_layer['material_type'] == 'TMD':
            # TMD-TMD interactions
            basis_size = len(lower_layer['unique_tags'])
            if basis_size >= 12:
                total_interactions += 64  # Orthorhombic
            else:
                # Hexagonal - classify and count
                from .input_generator import _classify_tmd_atoms, _should_interact_tmd
                lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                for tag_i, class_i, sym_i in lower_class:
                    for tag_j, class_j, sym_j in upper_class:
                        if _should_interact_tmd(class_i, class_j):
                            total_interactions += 1
        else:
            # hBN-hBN or hBN-TMD interactions
            lower_is_hbn = lower_layer['material_type'] == 'hBN'
            upper_is_hbn = upper_layer['material_type'] == 'hBN'
            
            if lower_is_hbn and upper_is_hbn:
                # hBN-hBN: all 4 interactions (B-B, B-N, N-B, N-N)
                total_interactions += 4
            elif lower_is_hbn and not upper_is_hbn:
                # hBN-TMD: both B and N interact with X_l (lower chalcogen) only
                # Classify TMD atoms to find X_l chalcogens
                from .input_generator import _classify_tmd_atoms
                upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                chalc_lower_count = sum(1 for tag, cls, sym in upper_class if cls == 'X_l')
                total_interactions += 2 * chalc_lower_count  # Both B and N from hBN
            elif not lower_is_hbn and upper_is_hbn:
                # TMD-hBN: X_u (upper chalcogen) interacts with both B and N
                # Classify TMD atoms to find X_u chalcogens
                from .input_generator import _classify_tmd_atoms
                lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                chalc_upper_count = sum(1 for tag, cls, sym in lower_class if cls == 'X_u')
                total_interactions += 2 * chalc_upper_count  # Both B and N from hBN
    
    print(f"  hBN layers: {n_hbn}, TMD layers: {n_tmd}")
    print(f"  Total interlayer interactions: {total_interactions}")
    
    # Calculate max atom type for proper array sizing
    max_atom_type = max(material_info['atom_types'])
    
    with open('lammps.in', 'w') as f:
        f.write("# General\n")
        f.write("units           metal\n")
        f.write("dimension       3\n")
        f.write("atom_style      atomic\n")
        f.write("neighbor        0.3 bin\n\n\n")
        
        f.write("# Structure\n")
        f.write("boundary        p p p\n")
        f.write(f"read_data       {structure_file}\n")
        for atype, mass in zip(material_info['atom_types'], material_info['masses']):
            f.write(f"mass {atype} {mass}\n")
        f.write("\n")
        
        # Pair style - separate sw/mod and tersoff for each layer
        f.write("# potential definitions\n")
        pair_style = "pair_style hybrid/overlay "
        
        # Add sw/mod for each TMD layer
        for info in material_info['layer_info']:
            if info['material_type'] == 'TMD':
                pair_style += "sw/mod "
        
        # Add tersoff for each hBN layer
        for info in material_info['layer_info']:
            if info['material_type'] == 'hBN':
                pair_style += "tersoff "
        
        # Add KC for all interlayer interactions
        pair_style += ' '.join([f'kolmogorov/crespi/z 14.0'] * total_interactions)
        pair_style += " lj/cut 10.0\n"
        f.write(pair_style)
        f.write(f"# Number of interlayer interactions = {total_interactions}\n\n")
        
        # Intralayer interactions
        f.write("# intralayer interactions\n")
        
        # TMD intralayer (SW/mod) - one instance per TMD layer
        tmd_layer_idx = 1
        for i, layer_info in enumerate(material_info['layer_info']):
            if layer_info['material_type'] == 'TMD':
                symbols_line = ['NULL'] * max_atom_type
                
                # Check if orthorhombic basis (12+ atom types)
                basis_size = len(layer_info['unique_tags'])
                is_ortho_basis = basis_size >= 12
                
                if is_ortho_basis:
                    # Orthorhombic: use numbered symbols (S1, W1, S2, etc.)
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
                f.write(f"pair_coeff * * sw/mod {tmd_layer_idx} {pot} {' '.join(symbols_line)}\n")
                tmd_layer_idx += 1
        
        # hBN intralayer (Tersoff) - one instance per hBN layer
        tersoff_layer_idx = 1
        for i, layer_info in enumerate(material_info['layer_info']):
            if layer_info['material_type'] == 'hBN':
                symbols_line = ['NULL'] * max_atom_type
                
                # Use unique_tags which already has correct types assigned
                # unique_tags = [type_B, type_N] from _analyze_system
                type_B = layer_info['unique_tags'][0]  # First is B
                type_N = layer_info['unique_tags'][1]  # Second is N
                
                symbols_line[type_B - 1] = 'B'
                symbols_line[type_N - 1] = 'N'
                
                pot = layer_info['intralayer_potential']
                f.write(f"pair_coeff * * tersoff {tersoff_layer_idx} {pot} {' '.join(symbols_line)}\n")
                tersoff_layer_idx += 1
        
        f.write("\npair_coeff * * lj/cut 0.0 3.4\n\n")
        
        # Interlayer interactions - all use KC
        f.write("# interlayer interaction\n")
        
        interaction_num = 1
        type_to_symbol = material_info['type_to_symbol']
        
        for i in range(n_layers - 1):
            lower_layer = material_info['layer_info'][i]
            upper_layer = material_info['layer_info'][i+1]
            
            # Determine potential file based on materials
            if lower_layer['material_type'] == 'TMD' and upper_layer['material_type'] == 'TMD':
                inter_pot = config.get('interlayer_potential', 'MoWSSe.KC')
            elif lower_layer['material_type'] == 'hBN' or upper_layer['material_type'] == 'hBN':
                inter_pot = 'SeBN.KC'  # hBN-TMD or hBN-hBN
            else:
                inter_pot = config.get('interlayer_potential', 'MoWSSe.KC')
            
            if lower_layer['material_type'] == 'TMD' and upper_layer['material_type'] == 'TMD':
                # TMD-TMD interactions
                basis_size = len(lower_layer['unique_tags'])
                if basis_size >= 12:
                    # Orthorhombic: 64 interactions
                    from .input_generator import _partition_chalcogens_ortho, _partition_metals_ortho
                    
                    lower_chalc_A, lower_chalc_B = _partition_chalcogens_ortho(lower_layer, type_to_symbol)
                    upper_chalc_A, upper_chalc_B = _partition_chalcogens_ortho(upper_layer, type_to_symbol)
                    lower_metals = _partition_metals_ortho(lower_layer, type_to_symbol)
                    upper_metals = _partition_metals_ortho(upper_layer, type_to_symbol)
                    
                    # Generate 64 interactions
                    for tag_i in lower_chalc_A:
                        for tag_j in upper_chalc_B:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = type_to_symbol[tag_i]
                            symbols_line[tag_j - 1] = type_to_symbol[tag_j]
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                    
                    for tag_i in lower_metals:
                        for tag_j in upper_metals:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = type_to_symbol[tag_i]
                            symbols_line[tag_j - 1] = type_to_symbol[tag_j]
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                    
                    for tag_i in lower_chalc_A:
                        for tag_j in upper_metals:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = type_to_symbol[tag_i]
                            symbols_line[tag_j - 1] = type_to_symbol[tag_j]
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                    
                    for tag_i in lower_metals:
                        for tag_j in upper_chalc_B:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = type_to_symbol[tag_i]
                            symbols_line[tag_j - 1] = type_to_symbol[tag_j]
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                else:
                    # Hexagonal TMD-TMD
                    from .input_generator import _classify_tmd_atoms, _should_interact_tmd
                    lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                    upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                    
                    for tag_i, class_i, sym_i in lower_class:
                        for tag_j, class_j, sym_j in upper_class:
                            if _should_interact_tmd(class_i, class_j):
                                symbols_line = ['NULL'] * max_atom_type
                                symbols_line[tag_i - 1] = sym_i
                                symbols_line[tag_j - 1] = sym_j
                                f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                                interaction_num += 1
            else:
                # hBN-hBN or hBN-TMD interactions
                TMD_CHALCOGENS = {'S', 'Se', 'Te'}
                
                # Get atom types for lower and upper layers
                lower_is_hbn = lower_layer['material_type'] == 'hBN'
                upper_is_hbn = upper_layer['material_type'] == 'hBN'
                
                if lower_is_hbn and not upper_is_hbn:
                    # hBN below TMD: only X_l (lower chalcogen) of TMD interacts with B and N
                    from .input_generator import _classify_tmd_atoms
                    upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                    
                    # Get hBN types and TMD X_l chalcogens
                    hbn_types = lower_layer['unique_tags']  # [B, N]
                    tmd_chalc_lower = [tag for tag, cls, sym in upper_class if cls == 'X_l']
                    
                    for tag_i in hbn_types:
                        for tag_j in tmd_chalc_lower:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = 'B'  # Always use B in symbols
                            symbols_line[tag_j - 1] = type_to_symbol[tag_j]
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                            
                elif upper_is_hbn and not lower_is_hbn:
                    # TMD below hBN: only X_u (upper chalcogen) of TMD interacts with B and N
                    from .input_generator import _classify_tmd_atoms
                    lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                    
                    # Get TMD X_u chalcogens and hBN types
                    tmd_chalc_upper = [tag for tag, cls, sym in lower_class if cls == 'X_u']
                    hbn_types = upper_layer['unique_tags']  # [B, N]
                    
                    for tag_i in tmd_chalc_upper:
                        for tag_j in hbn_types:
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = type_to_symbol[tag_i]
                            symbols_line[tag_j - 1] = 'B'  # Always use B in symbols
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
                            
                else:
                    # hBN-hBN: all interactions
                    lower_types = lower_layer['unique_tags']
                    upper_types = upper_layer['unique_tags']
                    
                    for tag_i in lower_types:
                        for tag_j in upper_types:
                            sym_i = type_to_symbol[tag_i]
                            sym_j = type_to_symbol[tag_j]
                            symbols_line = ['NULL'] * max_atom_type
                            symbols_line[tag_i - 1] = sym_i
                            symbols_line[tag_j - 1] = sym_j
                            f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                            interaction_num += 1
        
        f.write("\n# Optimize at 0 K\n")
        f.write("dump            1 all custom 400 dump.minimization id type x y z\n")
        f.write("thermo          500\n")
        f.write("min_style       fire\n")
        f.write("minimize        0.0 1.0e-4 1000000 1000000\n")
        f.write("write_dump all atom dump.Final\n")
        f.write('print "Done!"\n')


def generate_hbn_graphene_input(structure_file, material_info, config):
    """Generate lammps.in for mixed hBN + Graphene systems."""
    raise NotImplementedError("hBN + Graphene mixed systems not yet implemented")
