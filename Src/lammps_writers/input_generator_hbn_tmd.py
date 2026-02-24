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
                # hBN-hBN: all pairwise interactions between unique tags
                n_lower_hbn = len(lower_layer['unique_tags'])
                n_upper_hbn = len(upper_layer['unique_tags'])
                total_interactions += n_lower_hbn * n_upper_hbn
            elif lower_is_hbn and not upper_is_hbn:
                # hBN-TMD: all hBN tags interact with X_l (lower chalcogen) only
                from .input_generator import _classify_tmd_atoms
                upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                chalc_lower_count = sum(1 for tag, cls, sym in upper_class if cls == 'X_l')
                n_hbn_tags = len(lower_layer['unique_tags'])
                total_interactions += n_hbn_tags * chalc_lower_count
            elif not lower_is_hbn and upper_is_hbn:
                # TMD-hBN: X_u (upper chalcogen) interacts with all hBN tags
                from .input_generator import _classify_tmd_atoms
                lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                chalc_upper_count = sum(1 for tag, cls, sym in lower_class if cls == 'X_u')
                n_hbn_tags = len(upper_layer['unique_tags'])
                total_interactions += chalc_upper_count * n_hbn_tags
    
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
        
        type_to_symbol = material_info['type_to_symbol']
        
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
                
                # Map each unique tag to its symbol (B or N)
                for tag in layer_info['unique_tags']:
                    if tag in type_to_symbol:
                        symbols_line[tag - 1] = type_to_symbol[tag]
                    else:
                        # Fallback: first half are B, second half are N
                        n_tags = len(layer_info['unique_tags'])
                        idx_in_tags = layer_info['unique_tags'].index(tag)
                        symbols_line[tag - 1] = 'B' if idx_in_tags < n_tags // 2 else 'N'
                
                pot = layer_info['intralayer_potential']
                f.write(f"pair_coeff * * tersoff {tersoff_layer_idx} {pot} {' '.join(symbols_line)}\n")
                tersoff_layer_idx += 1
        
        f.write("\npair_coeff * * lj/cut 0.0 3.4\n\n")
        
        from .input_generator import (_classify_tmd_atoms, _should_interact_tmd,
                                        _partition_chalcogens_ortho, _partition_metals_ortho,
                                        _write_kc_pairs)
        
        interaction_num = 1
        
        # Collect layer pairs by interaction type
        tmd_tmd_pairs = []
        hbn_below_tmd_pairs = []  # hBN (lower) - TMD (upper)
        tmd_above_hbn_pairs = []  # TMD (lower) - hBN (upper)
        hbn_hbn_pairs = []
        
        for i in range(n_layers - 1):
            lower_layer = material_info['layer_info'][i]
            upper_layer = material_info['layer_info'][i+1]
            lower_is_hbn = lower_layer['material_type'] == 'hBN'
            upper_is_hbn = upper_layer['material_type'] == 'hBN'
            
            if not lower_is_hbn and not upper_is_hbn:
                tmd_tmd_pairs.append((i, i+1))
            elif lower_is_hbn and not upper_is_hbn:
                hbn_below_tmd_pairs.append((i, i+1))
            elif not lower_is_hbn and upper_is_hbn:
                tmd_above_hbn_pairs.append((i, i+1))
            else:
                hbn_hbn_pairs.append((i, i+1))
        
        # --- Section 1: TMD-TMD interlayer interactions ---
        if tmd_tmd_pairs:
            inter_pot = config.get('interlayer_potential', 'WSSe.KC')
            
            for li, ui in tmd_tmd_pairs:
                lower_layer = material_info['layer_info'][li]
                upper_layer = material_info['layer_info'][ui]
                basis_size = len(lower_layer['unique_tags'])
                
                if basis_size >= 12:
                    # Orthorhombic TMD-TMD
                    lower_chalc_A, lower_chalc_B = _partition_chalcogens_ortho(lower_layer, type_to_symbol, structure_file)
                    upper_chalc_A, upper_chalc_B = _partition_chalcogens_ortho(upper_layer, type_to_symbol, structure_file)
                    lower_metals = _partition_metals_ortho(lower_layer, type_to_symbol)
                    upper_metals = _partition_metals_ortho(upper_layer, type_to_symbol)
                    
                    lower_chalc_sym = type_to_symbol.get(lower_chalc_B[0], 'X') if lower_chalc_B else 'X'
                    upper_chalc_sym = type_to_symbol.get(upper_chalc_A[0], 'X') if upper_chalc_A else 'X'
                    metal_sym = type_to_symbol.get(lower_metals[0], 'M') if lower_metals else 'M'
                    
                    f.write(f"#------------------------------------------\n")
                    f.write(f"# {lower_chalc_sym}-{upper_chalc_sym} KC interactions.\n")
                    f.write(f"#------------------------------------------\n")
                    interaction_num = _write_kc_pairs(f, lower_chalc_B, upper_chalc_A, type_to_symbol, max_atom_type, inter_pot, interaction_num)
                    
                    f.write(f"\n#------------------------------------------\n")
                    f.write(f"# {metal_sym}-{metal_sym} KC interactions.\n")
                    f.write(f"#------------------------------------------\n")
                    interaction_num = _write_kc_pairs(f, lower_metals, upper_metals, type_to_symbol, max_atom_type, inter_pot, interaction_num)
                    
                    f.write(f"\n#------------------------------------------\n")
                    f.write(f"# {lower_chalc_sym}-{metal_sym} KC interactions.\n")
                    f.write(f"#------------------------------------------\n")
                    interaction_num = _write_kc_pairs(f, lower_chalc_B, upper_metals, type_to_symbol, max_atom_type, inter_pot, interaction_num)
                    
                    f.write(f"\n#------------------------------------------\n")
                    f.write(f"# {metal_sym}-{upper_chalc_sym} KC interactions.\n")
                    f.write(f"#------------------------------------------\n")
                    interaction_num = _write_kc_pairs(f, lower_metals, upper_chalc_A, type_to_symbol, max_atom_type, inter_pot, interaction_num)
                else:
                    # Hexagonal TMD-TMD
                    lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                    upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                    
                    f.write(f"#------------------------------------------\n")
                    f.write(f"# TMD-TMD KC interactions.\n")
                    f.write(f"#------------------------------------------\n")
                    for tag_i, class_i, sym_i in lower_class:
                        for tag_j, class_j, sym_j in upper_class:
                            if _should_interact_tmd(class_i, class_j):
                                symbols_line = ['NULL'] * max_atom_type
                                symbols_line[tag_i - 1] = sym_i
                                symbols_line[tag_j - 1] = sym_j
                                f.write(f"pair_coeff {tag_i} {tag_j} kolmogorov/crespi/z {interaction_num} {inter_pot} {' '.join(symbols_line)}\n")
                                interaction_num += 1
        
        # --- Section 2: hBN (bottom) - TMD interlayer interactions ---
        if hbn_below_tmd_pairs:
            for li, ui in hbn_below_tmd_pairs:
                lower_layer = material_info['layer_info'][li]  # hBN
                upper_layer = material_info['layer_info'][ui]  # TMD
                
                upper_class = _classify_tmd_atoms(upper_layer, structure_file)
                tmd_chalc_lower = [tag for tag, cls, sym in upper_class if cls == 'X_l']
                
                chalc_sym = type_to_symbol.get(tmd_chalc_lower[0], 'Se') if tmd_chalc_lower else 'Se'
                inter_pot = f"{chalc_sym}BN.KC"
                hbn_override = {tag: 'B' for tag in lower_layer['unique_tags']}
                
                f.write(f"\n#------------------------------------------\n")
                f.write(f"# {chalc_sym}-BN KC interactions (hBN bottom layer).\n")
                f.write(f"#------------------------------------------\n")
                interaction_num = _write_kc_pairs(f, lower_layer['unique_tags'], tmd_chalc_lower,
                                                 type_to_symbol, max_atom_type, inter_pot, interaction_num,
                                                 sym_override=hbn_override)
        
        # --- Section 3: TMD - hBN (top) interlayer interactions ---
        if tmd_above_hbn_pairs:
            for li, ui in tmd_above_hbn_pairs:
                lower_layer = material_info['layer_info'][li]  # TMD
                upper_layer = material_info['layer_info'][ui]  # hBN
                
                lower_class = _classify_tmd_atoms(lower_layer, structure_file)
                tmd_chalc_upper = [tag for tag, cls, sym in lower_class if cls == 'X_u']
                
                chalc_sym = type_to_symbol.get(tmd_chalc_upper[0], 'Se') if tmd_chalc_upper else 'Se'
                inter_pot = f"{chalc_sym}BN.KC"
                hbn_override = {tag: 'B' for tag in upper_layer['unique_tags']}
                
                f.write(f"\n#------------------------------------------\n")
                f.write(f"# {chalc_sym}-BN KC interactions (hBN top layer).\n")
                f.write(f"#------------------------------------------\n")
                interaction_num = _write_kc_pairs(f, tmd_chalc_upper, upper_layer['unique_tags'],
                                                 type_to_symbol, max_atom_type, inter_pot, interaction_num,
                                                 sym_override=hbn_override)
        
        # --- Section 4: hBN-hBN interlayer interactions (if any) ---
        if hbn_hbn_pairs:
            f.write(f"\n#------------------------------------------\n")
            f.write(f"# BN-BN KC interactions.\n")
            f.write(f"#------------------------------------------\n")
            for li, ui in hbn_hbn_pairs:
                lower_layer = material_info['layer_info'][li]
                upper_layer = material_info['layer_info'][ui]
                inter_pot = 'BNBN.KC'
                interaction_num = _write_kc_pairs(f, lower_layer['unique_tags'], upper_layer['unique_tags'],
                                                 type_to_symbol, max_atom_type, inter_pot, interaction_num)
        
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
