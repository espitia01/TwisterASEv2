import numpy as np
import sys

def parse_input_file(filename):
    """
    Parse main twisterase.inp file.
    Returns dict with configuration parameters.
    """
    config = {
        'n_layers': None,
        'hex_lattice': True,
        'i_value': None,
        'mn_values': None,
        'lattice_parameters': None,
        'superlattice_vectors': None,
        'write_lammps': False,
        'interlayer_potential': None
    }
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"{filename} not found!")
        sys.exit(1)
    
    for iline, line in enumerate(lines):
        keyword, _, value = line.partition('#')[0].partition('=')
        keyword = keyword.strip()
        
        if "n_layers" in keyword:
            config['n_layers'] = eval(value)
        elif "hex_lattice" in keyword:
            config['hex_lattice'] = eval(value)
        elif "i_value" in keyword:
            config['i_value'] = eval(value)
        elif "mn_values" in keyword:
            config['mn_values'] = eval(value)
        elif "lattice_parameters" in keyword:
            config['lattice_parameters'] = np.array(eval(value))
        elif "superlattice_vectors_block" in keyword:
            sl_vectors = np.zeros((3, 3))
            for iv in range(3):
                sl_vectors[iv] = [float(x) for x in lines[iline+1+iv].split()]
            config['superlattice_vectors'] = sl_vectors
        elif "write_lammps" in keyword:
            config['write_lammps'] = eval(value)
        elif "interlayer_potential" in keyword:
            config['interlayer_potential'] = eval(value)
    
    if config['n_layers'] is None:
        print("n_layers not set in input file!")
        sys.exit(1)
    
    return config


def parse_layer_file(filename):
    """
    Parse individual layer input file.
    Returns dict with layer parameters.
    """
    layer_data = {
        'twist_angle': 0.0,
        'lattice_parameters': None,
        'lattice_vectors': None,
        'atom_positions': None,
        'atom_symbols': None,
        'atom_tags': None,
        'orthocell_tm': None,
        'atom_positions_ortho': None,
        'atom_symbols_ortho': None,
        'atom_tags_ortho': None,
        'translate_z': 0.0,
        'intralayer_potential': None
    }
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"{filename} not found!")
        sys.exit(1)
    
    for iline, line in enumerate(lines):
        keyword, _, value = line.partition('#')[0].partition('=')
        keyword = keyword.strip()
        
        if "twist_angle" in keyword:
            layer_data['twist_angle'] = eval(value)
        elif "lattice_parameters" in keyword:
            layer_data['lattice_parameters'] = np.array(eval(value))
        elif "translate_z" in keyword:
            layer_data['translate_z'] = eval(value)
        elif "intralayer_potential" in keyword:
            layer_data['intralayer_potential'] = eval(value)
        elif "lattice_vectors_block" in keyword:
            lv = np.zeros((3, 3))
            for iv in range(3):
                lv[iv] = [float(x) for x in lines[iline+1+iv].split()]
            layer_data['lattice_vectors'] = lv
        elif "orthocell_transformation_matrix" in keyword:
            tm = np.zeros((3, 3))
            for iv in range(3):
                tm[iv] = [float(x) for x in lines[iline+1+iv].split()]
            layer_data['orthocell_tm'] = tm
        elif "start_atom_positions_block" in keyword:
            positions = []
            symbols = []
            tags = []
            for at_line in lines[iline+1:]:
                if "end_atom_positions_block" in at_line:
                    break
                parts = at_line.split()
                symbols.append(parts[0])
                positions.append([float(x) for x in parts[1:4]])
                if len(parts) > 4:
                    tags.append(int(parts[4]))
            layer_data['atom_positions'] = np.array(positions)
            layer_data['atom_symbols'] = symbols
            layer_data['atom_tags'] = tags if tags else None
        elif "start_atom_positions_ortho_block" in keyword:
            positions = []
            symbols = []
            tags = []
            for at_line in lines[iline+1:]:
                if "end_atom_positions_ortho_block" in at_line:
                    break
                parts = at_line.split()
                symbols.append(parts[0])
                positions.append([float(x) for x in parts[1:4]])
                if len(parts) > 4:
                    tags.append(int(parts[4]))
            layer_data['atom_positions_ortho'] = np.array(positions)
            layer_data['atom_symbols_ortho'] = symbols
            layer_data['atom_tags_ortho'] = tags if tags else None
    
    return layer_data
