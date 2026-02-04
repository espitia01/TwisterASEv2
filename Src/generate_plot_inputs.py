#!/usr/bin/env python3
"""
generate_plot_inputs.py - Generate input files for plotting scripts

Reads cutpos.inp and layer*.inp files to generate 'input' files
in each analysis subdirectory with proper lattice vectors and parameters.
"""

import os
import re
import numpy as np


def get_superlattice_vectors_from_structure():
    """Get actual superlattice vectors from relaxed_structure.cif."""
    # Read CIF file manually to avoid scipy import issues
    import warnings
    warnings.filterwarnings('ignore')
    
    try:
        from ase.io import read
        # Try to read from relaxed_structure.cif first
        if os.path.exists('relaxed_structure.cif'):
            atoms = read('relaxed_structure.cif')
            cell = atoms.get_cell()
            # Return cell vectors in Angstroms
            return cell[:3]
        elif os.path.exists('superlattice.cif'):
            atoms = read('superlattice.cif')
            cell = atoms.get_cell()
            return cell[:3]
        else:
            raise FileNotFoundError("Could not find relaxed_structure.cif or superlattice.cif")
    except (ValueError, ImportError) as e:
        # Fallback: parse CIF manually
        print(f"Warning: ASE import issue ({e}), using manual CIF parsing")
        return parse_cif_cell_manually()


def parse_cif_cell_manually():
    """Manually parse CIF file to extract cell parameters."""
    cif_file = 'relaxed_structure.cif' if os.path.exists('relaxed_structure.cif') else 'superlattice.cif'
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    # Extract cell parameters from CIF
    cell_length_a = cell_length_b = cell_length_c = None
    cell_angle_alpha = cell_angle_beta = cell_angle_gamma = None
    
    for line in lines:
        if '_cell_length_a' in line:
            cell_length_a = float(line.split()[1])
        elif '_cell_length_b' in line:
            cell_length_b = float(line.split()[1])
        elif '_cell_length_c' in line:
            cell_length_c = float(line.split()[1])
        elif '_cell_angle_alpha' in line:
            cell_angle_alpha = float(line.split()[1])
        elif '_cell_angle_beta' in line:
            cell_angle_beta = float(line.split()[1])
        elif '_cell_angle_gamma' in line:
            cell_angle_gamma = float(line.split()[1])
    
    # Convert to Cartesian vectors
    a = cell_length_a
    b = cell_length_b
    c = cell_length_c
    alpha = np.radians(cell_angle_alpha)
    beta = np.radians(cell_angle_beta)
    gamma = np.radians(cell_angle_gamma)
    
    # Standard conversion
    A1 = np.array([a, 0, 0])
    A2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    A3 = np.array([cx, cy, cz])
    
    return np.array([A1, A2, A3])


def parse_layer_inp(filepath):
    """Parses a layer#.inp file to extract alat and atom of interest."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract lattice_parameters (alat)
    alat_match = re.search(r'lattice_parameters\s*=\s*\[(.*?)\]', content, re.IGNORECASE)
    if not alat_match:
        raise ValueError(f"Could not find 'lattice_parameters' in {filepath}")
    alat = float(alat_match.group(1).split(',')[0])
    
    # Extract atom of interest (TM, C, or B)
    atom_block_match = re.search(r'start_atom_positions_block\s*([\s\S]*?)end_atom_positions_block', content, re.IGNORECASE)
    if not atom_block_match:
        raise ValueError(f"Could not find atom positions block in {filepath}")
    
    atom_lines = atom_block_match.group(1).strip().split('\n')
    atom_of_interest = None
    for line in atom_lines:
        parts = line.split()
        if len(parts) > 0:
            atom_symbol = parts[0]
            if atom_symbol in ['W', 'Mo', 'C', 'B', 'N', 'Bi', 'Se', 'S', 'Te']:
                atom_of_interest = atom_symbol
                break
    
    if not atom_of_interest:
        raise ValueError(f"Could not determine atom of interest from {filepath}")
    
    return alat, atom_of_interest


def write_plot_input(f, alat, A1, A2, A3, spcs, repeat_units=(2, 2), shift_origin=(2, 2)):
    """Write a plot input file with dynamically calculated limits."""
    f.write(f"alat = {alat}\n")
    f.write(f"A1 = {A1[0]:.14f} {A1[1]:.14f} {A1[2]:.14f}\n")
    f.write(f"A2 = {A2[0]:.14f} {A2[1]:.14f} {A2[2]:.14f}\n")
    f.write(f"A3 = {A3[0]:.14f} {A3[1]:.14f} {A3[2]:.14f}\n")
    f.write(f"compute_spcs = [{spcs[0]} {spcs[1]}]\n")
    f.write("rotate_plot = False\n")
    f.write("rotate_angle = 0\n")
    
    # Calculate xlim and ylim based on actual repeat_units and lattice vectors
    # The data positions are: pos = i*A1 + j*A2 - origin
    # where origin = shift_origin[0]*A1 + shift_origin[1]*A2
    # With padding=2, i ranges from -padding to repeat_units[0] + padding
    padding = 2
    repeat_x, repeat_y = repeat_units
    shift_x, shift_y = shift_origin
    
    # Calculate the origin shift
    origin = shift_x * A1 * alat + shift_y * A2 * alat
    
    # Calculate the extent of the supercell with padding
    min_i = -padding
    max_i = repeat_x + padding
    min_j = -padding
    max_j = repeat_y + padding
    
    # Calculate corner positions in Cartesian coordinates
    # pos = i*A1 + j*A2 - origin
    corners_x = []
    corners_y = []
    for i in [min_i, max_i]:
        for j in [min_j, max_j]:
            pos = i * A1 * alat + j * A2 * alat - origin
            corners_x.append(pos[0])
            corners_y.append(pos[1])
    
    # Get min/max with small padding
    min_x = min(corners_x)
    max_x = max(corners_x)
    min_y = min(corners_y)
    max_y = max(corners_y)
    
    # Add 5% padding for visualization
    x_range = max_x - min_x
    y_range = max_y - min_y
    xlim_min = int(min_x - 0.05 * x_range)
    xlim_max = int(max_x + 0.05 * x_range)
    ylim_min = int(min_y - 0.05 * y_range)
    ylim_max = int(max_y + 0.05 * y_range)
    
    f.write(f"xlim = {xlim_min} {xlim_max}\n")
    f.write(f"ylim = {ylim_min} {ylim_max}\n")
    f.write(f"repeat_units = {repeat_units[0]} {repeat_units[1]}\n")
    f.write("place_dots = False\n")
    f.write("fontsize = 20\n")
    f.write("tick_labelsize = 16\n")
    f.write(f"shift_origin = {shift_origin[0]} {shift_origin[1]}\n")


def main():
    """Main function to generate all input files."""
    print("Generating plot input files...")
    base_path = '.'
    
    # Get actual superlattice vectors from the relaxed structure
    superlattice_vectors = get_superlattice_vectors_from_structure()
    A1, A2, A3 = superlattice_vectors
    
    # Get layer info
    layer_files = sorted([f for f in os.listdir(base_path) if re.match(r'layer\d+\.inp', f)])
    n_layers = len(layer_files)
    layer_params = [parse_layer_inp(os.path.join(base_path, f)) for f in layer_files]
    
    # --- Generate Interlayer Spacing Inputs (only for multi-layer systems) ---
    if n_layers > 1:
        ils_path = os.path.join(base_path, 'InterlayerSpacingMap')
        for i in range(n_layers - 1):
            # For ILS, use cut vectors directly with alat = 1.0
            spcs = [layer_params[i][1], layer_params[i+1][1]]
            
            input_dir = os.path.join(ils_path, f'Layer_{i+1}-{i+2}')
            with open(os.path.join(input_dir, 'input'), 'w') as f:
                write_plot_input(f, 1.0, A1, A2, A3, spcs, repeat_units=(2, 2), shift_origin=(0, 0))
            print(f"  - Wrote input for {input_dir}")
    else:
        print("  - Skipping interlayer spacing inputs for monolayer system")
    
    # --- Generate Strain Inputs ---
    strain_path = os.path.join(base_path, 'StrainMap')
    for i in range(n_layers):
        alat = layer_params[i][0]
        atom_interest = layer_params[i][1]
        # For strain, divide cut vectors by the layer's lattice constant
        scaled_A1 = A1 / alat
        scaled_A2 = A2 / alat
        scaled_A3 = A3 / alat
        spcs = [atom_interest, atom_interest]
        
        input_dir = os.path.join(strain_path, f'Layer_{i+1}')
        with open(os.path.join(input_dir, 'input'), 'w') as f:
            write_plot_input(f, alat, scaled_A1, scaled_A2, scaled_A3, spcs, repeat_units=(2, 2), shift_origin=(1, 1))
        print(f"  - Wrote input for {input_dir}")
    
    print("\nInput file generation complete.")


if __name__ == "__main__":
    main()
