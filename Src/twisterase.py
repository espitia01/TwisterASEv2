import numpy as np
from file_io import parse_input_file
from layer import Layer
from transformations import angle_tm_using_i, angle_tm_using_mn, compute_ortho_supercell_matrix
from ase.io import write


def main():
    config = parse_input_file('twisterase.inp')
    
    n_layers = config['n_layers']
    sl_vectors = config['superlattice_vectors']
    sl_alats = config['lattice_parameters']
    i_value = config['i_value']
    mn_values = config['mn_values']
    hex_lattice = config['hex_lattice']
    
    if sl_alats is not None and sl_vectors is not None:
        sl_vectors_cart = sl_vectors * sl_alats
    else:
        sl_vectors_cart = None
    
    twist_angle_from_i = None
    if hex_lattice and i_value is not None:
        layer1 = Layer('layer1.inp')
        twist_angle_from_i, sc_matrix = angle_tm_using_i(i_value)
        sl_vectors_cart = np.dot(sc_matrix, layer1.unitcell.cell)
        print(f"Using i_value={i_value}, twist angle={twist_angle_from_i:.4f} degrees")
    elif hex_lattice and mn_values is not None:
        layer1 = Layer('layer1.inp')
        twist_angle_from_i, sc_matrix = angle_tm_using_mn(mn_values[0], mn_values[1])
        sl_vectors_cart = np.dot(sc_matrix, layer1.unitcell.cell)
        print(f"Using mn_values={mn_values}, twist angle={twist_angle_from_i:.4f} degrees")
    
    layers = []
    has_orthocell = False
    
    for ilayer in range(n_layers):
        layer_file = f'layer{ilayer+1}.inp'
        layer = Layer(layer_file)
        
        # Override twist angles if i_value or mn_values is specified
        if i_value is not None or mn_values is not None:
            if ilayer == 0:
                layer.set_twist_angle(0.0)
                param_name = "i_value" if i_value is not None else "mn_values"
                print(f"Layer {ilayer+1}: twist angle set to 0.0 degrees (from {param_name})")
            elif ilayer == 1:
                layer.set_twist_angle(twist_angle_from_i)
                param_name = "i_value" if i_value is not None else "mn_values"
                print(f"Layer {ilayer+1}: twist angle set to {twist_angle_from_i:.4f} degrees (from {param_name})")
        
        layer.create_supercell(sl_vectors_cart)
        layers.append(layer)
        
        if layer.has_orthocell:
            has_orthocell = True
    
    superlattice = layers[0].supercell.copy()
    for ilayer in range(1, n_layers):
        superlattice += layers[ilayer].supercell
    
    write('superlattice.cif', superlattice)
    print(f"Written superlattice.cif with {len(superlattice)} atoms")
    
    if has_orthocell:
        hex_sl_vectors = superlattice.get_cell()
        ortho_unitcell_vectors = layers[0].unitcell_ortho.cell
        sl_vectors_cart_ortho = compute_ortho_supercell_matrix(
            hex_sl_vectors, ortho_unitcell_vectors
        )
        
        for layer in layers:
            layer.create_supercell_ortho(sl_vectors_cart_ortho)
        
        superlattice_ortho = layers[0].supercell_ortho.copy()
        for ilayer in range(1, n_layers):
            superlattice_ortho += layers[ilayer].supercell_ortho
        
        #write('superlattice_ortho.cif', superlattice_ortho)
        #print(f"Written superlattice_ortho.cif with {len(superlattice_ortho)} atoms")
        
        if config['write_lammps']:
            from lammps_writers import write_structure_file, generate_lammps_input
            material_info = write_structure_file('structure_ortho.lammps', layers, use_orthocell=True)
            print(f"Written structure_ortho.lammps (Materials: {', '.join(material_info['layer_types'])})")
            generate_lammps_input('structure_ortho.lammps', layers, config, use_orthocell=True)
    elif config['write_lammps']:
        from lammps_writers import write_structure_file, generate_lammps_input
        material_info = write_structure_file('structure.lammps', layers, use_orthocell=False)
        print(f"Written structure.lammps (Materials: {', '.join(material_info['layer_types'])})")
        generate_lammps_input('structure.lammps', layers, config, use_orthocell=False)


if __name__ == '__main__':
    main()
