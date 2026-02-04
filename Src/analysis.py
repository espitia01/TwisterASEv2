#!/usr/bin/env python3
"""
analysis.py - Analyze interlayer spacing and strain in relaxed structures

Reads relaxed structures and layer CIF files to compute:
- Interlayer spacing distribution
- In-plane strain distribution
- Visualizations
"""

import numpy as np
import sys
from ase.io import read
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree


def read_layer_coords(filename):
    """
    Read layer coordinates from .dat file (faster than CIF).
    
    Args:
        filename: Path to layer_#_coords.dat file
    
    Returns:
        numpy array of positions (N x 3)
    """
    positions = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                # Format: symbol x y z
                pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                positions.append(pos)
    return np.array(positions)


def compute_interlayer_spacing(relaxed_file, layer_files):
    """
    Compute interlayer spacing between layers.
    
    Args:
        relaxed_file: Path to relaxed structure (CIF or dump.Final)
        layer_files: List of paths to layer files (CIF or .dat coords)
    
    Returns:
        dict with spacing statistics
    """
    print(f"\n{'='*60}")
    print("Computing Interlayer Spacing")
    print(f"{'='*60}\n")
    
    # Read relaxed structure
    if relaxed_file.endswith('.cif'):
        atoms = read(relaxed_file)
        positions = atoms.get_positions()
    else:
        # Assume dump.Final format
        atoms = read(relaxed_file, format='lammps-dump-text', index=-1)
        positions = atoms.get_positions()
    
    # Read layer files (support both CIF and .dat formats)
    layer_positions = []
    for layer_file in layer_files:
        if layer_file.endswith('.dat'):
            # Fast: read .dat coordinate file directly
            pos = read_layer_coords(layer_file)
            layer_positions.append(pos)
            print(f"  Read {layer_file}: {len(pos)} atoms (from .dat)")
        else:
            # Slower: read CIF file
            layer = read(layer_file)
            layer_positions.append(layer.get_positions())
            print(f"  Read {layer_file}: {len(layer)} atoms (from CIF)")
    
    n_layers = len(layer_positions)
    print(f"\nAnalyzing {n_layers} layers")
    print(f"Total atoms in relaxed structure: {len(positions)}")
    
    # Compute average z-position for each layer
    layer_z_avg = []
    for i, pos in enumerate(layer_positions):
        avg_z = np.mean(pos[:, 2])
        layer_z_avg.append(avg_z)
        print(f"  Layer {i+1}: <z> = {avg_z:.4f} Å ({len(pos)} atoms)")
    
    # Compute interlayer spacings
    spacings = []
    for i in range(n_layers - 1):
        spacing = layer_z_avg[i+1] - layer_z_avg[i]
        spacings.append(spacing)
        print(f"\nInterlayer spacing {i+1}→{i+2}: {spacing:.4f} Å")
    
    # Compute detailed spacing distribution
    if n_layers == 2:
        layer1_pos = layer_positions[0]
        layer2_pos = layer_positions[1]
        
        # For each atom in layer 1, find nearest neighbor in layer 2
        tree = cKDTree(layer2_pos[:, :2])  # Only x,y coordinates
        local_spacings = []
        
        for pos1 in layer1_pos:
            dist, idx = tree.query(pos1[:2])
            pos2 = layer2_pos[idx]
            local_spacing = pos2[2] - pos1[2]
            local_spacings.append(local_spacing)
        
        local_spacings = np.array(local_spacings)
        
        print(f"\nLocal spacing statistics:")
        print(f"  Mean:   {np.mean(local_spacings):.4f} Å")
        print(f"  Std:    {np.std(local_spacings):.4f} Å")
        print(f"  Min:    {np.min(local_spacings):.4f} Å")
        print(f"  Max:    {np.max(local_spacings):.4f} Å")
        print(f"  Range:  {np.max(local_spacings) - np.min(local_spacings):.4f} Å")
        
        return {
            'average_spacing': spacings[0],
            'local_spacings': local_spacings,
            'layer_z_positions': layer_z_avg,
            'mean': np.mean(local_spacings),
            'std': np.std(local_spacings),
            'min': np.min(local_spacings),
            'max': np.max(local_spacings)
        }
    
    return {
        'spacings': spacings,
        'layer_z_positions': layer_z_avg
    }


def compute_strain(initial_cif, relaxed_cif):
    """
    Compute in-plane strain by comparing initial and relaxed structures.
    
    Args:
        initial_cif: Path to initial structure CIF
        relaxed_cif: Path to relaxed structure CIF
    
    Returns:
        dict with strain statistics
    """
    print(f"\n{'='*60}")
    print("Computing In-Plane Strain")
    print(f"{'='*60}\n")
    
    initial = read(initial_cif)
    relaxed = read(relaxed_cif)
    
    initial_cell = initial.get_cell()
    relaxed_cell = relaxed.get_cell()
    
    # Compute lattice parameter changes
    a_initial = np.linalg.norm(initial_cell[0])
    b_initial = np.linalg.norm(initial_cell[1])
    
    a_relaxed = np.linalg.norm(relaxed_cell[0])
    b_relaxed = np.linalg.norm(relaxed_cell[1])
    
    strain_a = (a_relaxed - a_initial) / a_initial * 100
    strain_b = (b_relaxed - b_initial) / b_initial * 100
    
    print(f"Lattice parameter a:")
    print(f"  Initial:  {a_initial:.6f} Å")
    print(f"  Relaxed:  {a_relaxed:.6f} Å")
    print(f"  Strain:   {strain_a:.4f} %")
    
    print(f"\nLattice parameter b:")
    print(f"  Initial:  {b_initial:.6f} Å")
    print(f"  Relaxed:  {b_relaxed:.6f} Å")
    print(f"  Strain:   {strain_b:.4f} %")
    
    # Compute local strain (bond length changes)
    if len(initial) == len(relaxed):
        initial_pos = initial.get_positions()
        relaxed_pos = relaxed.get_positions()
        
        # Compute nearest neighbor distances
        tree_initial = cKDTree(initial_pos[:, :2])
        tree_relaxed = cKDTree(relaxed_pos[:, :2])
        
        local_strains = []
        for i in range(len(initial)):
            # Find 3 nearest neighbors (excluding self)
            dist_init, idx_init = tree_initial.query(initial_pos[i, :2], k=4)
            dist_relax, idx_relax = tree_relaxed.query(relaxed_pos[i, :2], k=4)
            
            # Average of nearest neighbor distances (skip first which is self)
            avg_dist_init = np.mean(dist_init[1:])
            avg_dist_relax = np.mean(dist_relax[1:])
            
            if avg_dist_init > 0:
                local_strain = (avg_dist_relax - avg_dist_init) / avg_dist_init * 100
                local_strains.append(local_strain)
        
        local_strains = np.array(local_strains)
        
        print(f"\nLocal strain statistics:")
        print(f"  Mean:   {np.mean(local_strains):.4f} %")
        print(f"  Std:    {np.std(local_strains):.4f} %")
        print(f"  Min:    {np.min(local_strains):.4f} %")
        print(f"  Max:    {np.max(local_strains):.4f} %")
        
        return {
            'strain_a': strain_a,
            'strain_b': strain_b,
            'local_strains': local_strains,
            'mean': np.mean(local_strains),
            'std': np.std(local_strains)
        }
    
    return {
        'strain_a': strain_a,
        'strain_b': strain_b
    }


def plot_interlayer_spacing(spacing_data, output_file='interlayer_spacing.png'):
    """Plot interlayer spacing distribution."""
    if 'local_spacings' not in spacing_data:
        print("No local spacing data available for plotting")
        return
    
    local_spacings = spacing_data['local_spacings']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1.hist(local_spacings, bins=30, edgecolor='black', alpha=0.7)
    ax1.axvline(spacing_data['mean'], color='red', linestyle='--', 
                label=f"Mean: {spacing_data['mean']:.4f} Å")
    ax1.set_xlabel('Interlayer Spacing (Å)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Interlayer Spacing Distribution', fontsize=14)
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Box plot
    ax2.boxplot(local_spacings, vert=True)
    ax2.set_ylabel('Interlayer Spacing (Å)', fontsize=12)
    ax2.set_title('Interlayer Spacing Statistics', fontsize=14)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Saved plot: {output_file}")
    plt.close()


def plot_strain(strain_data, output_file='strain_distribution.png'):
    """Plot strain distribution."""
    if 'local_strains' not in strain_data:
        print("No local strain data available for plotting")
        return
    
    local_strains = strain_data['local_strains']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1.hist(local_strains, bins=30, edgecolor='black', alpha=0.7, color='orange')
    ax1.axvline(strain_data['mean'], color='red', linestyle='--',
                label=f"Mean: {strain_data['mean']:.4f} %")
    ax1.set_xlabel('Local Strain (%)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('In-Plane Strain Distribution', fontsize=14)
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Box plot
    ax2.boxplot(local_strains, vert=True)
    ax2.set_ylabel('Local Strain (%)', fontsize=12)
    ax2.set_title('Strain Statistics', fontsize=14)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ Saved plot: {output_file}")
    plt.close()


def main():
    """Main analysis function."""
    if len(sys.argv) < 2:
        print("Usage: python analysis.py <command> [options]")
        print("\nCommands:")
        print("  spacing <relaxed> <layer1> <layer2> ...")
        print("    - Compute interlayer spacing")
        print("\n  strain <initial> <relaxed>")
        print("    - Compute in-plane strain")
        print("\n  all <initial> <relaxed> <layer1> <layer2> ...")
        print("    - Run all analyses")
        print("\nFile formats supported:")
        print("  - CIF files: relaxed_structure.cif, layer_1.cif, etc.")
        print("  - LAMMPS dump: dump.Final")
        print("  - Coordinate files: layer_1_coords.dat, layer_2_coords.dat (FAST!)")
        print("\nExamples:")
        print("  # Using .dat files (fastest):")
        print("  python analysis.py spacing dump.Final layer_1_coords.dat layer_2_coords.dat")
        print("\n  # Using CIF files:")
        print("  python analysis.py spacing relaxed_structure.cif layer_1.cif layer_2.cif")
        print("\n  # Strain analysis:")
        print("  python analysis.py strain superlattice.cif relaxed_structure.cif")
        print("\n  # Complete analysis with .dat files:")
        print("  python analysis.py all superlattice.cif dump.Final layer_1_coords.dat layer_2_coords.dat")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == 'spacing':
        if len(sys.argv) < 4:
            print("❌ Error: Need relaxed CIF and at least 2 layer CIF files")
            sys.exit(1)
        
        relaxed_cif = sys.argv[2]
        layer_cifs = sys.argv[3:]
        
        spacing_data = compute_interlayer_spacing(relaxed_cif, layer_cifs)
        plot_interlayer_spacing(spacing_data)
        
    elif command == 'strain':
        if len(sys.argv) < 4:
            print("❌ Error: Need initial and relaxed CIF files")
            sys.exit(1)
        
        initial_cif = sys.argv[2]
        relaxed_cif = sys.argv[3]
        
        strain_data = compute_strain(initial_cif, relaxed_cif)
        plot_strain(strain_data)
        
    elif command == 'all':
        if len(sys.argv) < 5:
            print("❌ Error: Need initial CIF, relaxed CIF, and layer CIF files")
            sys.exit(1)
        
        initial_cif = sys.argv[2]
        relaxed_cif = sys.argv[3]
        layer_cifs = sys.argv[4:]
        
        spacing_data = compute_interlayer_spacing(relaxed_cif, layer_cifs)
        plot_interlayer_spacing(spacing_data)
        
        strain_data = compute_strain(initial_cif, relaxed_cif)
        plot_strain(strain_data)
        
        print(f"\n{'='*60}")
        print("Analysis Complete!")
        print(f"{'='*60}")
        
    else:
        print(f"❌ Unknown command: {command}")
        sys.exit(1)


if __name__ == '__main__':
    main()
