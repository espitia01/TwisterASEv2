#!/usr/bin/env python3
"""
setup_analysis_dirs.py - Setup directory structure for analysis

Creates organized folders:
- InterlayerSpacingMap/Layer_X-Y/ (for bilayers and multilayers)
- StrainMap/Layer_X/ (for all layers)

Copies layer coordinate files to appropriate directories.
"""

import os
import shutil
import glob


def setup_analysis_dirs(base_path='.'):
    """
    Sets up the directory structure for interlayer spacing and strain analysis.
    
    Creates 'InterlayerSpacingMap' and 'StrainMap' directories,
    then populates them with subdirectories and atom position files.
    
    Args:
        base_path (str): The root directory where the analysis folders will be created.
    """
    print("Setting up analysis directories...")
    
    # Find all layer position files
    pos_files = sorted(glob.glob(os.path.join(base_path, 'layer_*_coords.dat')))
    n_layers = len(pos_files)
    
    if n_layers == 0:
        print("Error: No 'layer_*_coords.dat' files found in the current directory.")
        print("Please run cutpos.py first to extract layer coordinates.")
        return
    
    print(f"Found {n_layers} layer position files.")
    
    # --- Create Interlayer Spacing Directories (only for multi-layer systems) ---
    if n_layers > 1:
        ils_path = os.path.join(base_path, 'InterlayerSpacingMap')
        os.makedirs(ils_path, exist_ok=True)
        print(f"Created directory: {ils_path}")
        
        print(f"Setting up Interlayer Spacing directory: {ils_path}")
        
        for i in range(n_layers - 1):
            layer_dir = os.path.join(ils_path, f'Layer_{i+1}-{i+2}')
            os.makedirs(layer_dir, exist_ok=True)
            
            # Copy and rename position files
            # pos_l = lower layer, pos_u = upper layer
            shutil.copy(pos_files[i], os.path.join(layer_dir, 'pos_l'))
            shutil.copy(pos_files[i+1], os.path.join(layer_dir, 'pos_u'))
            print(f"  - Created subdirectory {layer_dir} and copied position files.")
    else:
        print("Skipping interlayer spacing setup for monolayer system.")
    
    # --- Create Strain Directories ---
    strain_path = os.path.join(base_path, 'StrainMap')
    os.makedirs(strain_path, exist_ok=True)
    print(f"Created directory: {strain_path}")
    
    print(f"Setting up Strain Map directory: {strain_path}")
    
    for i in range(n_layers):
        layer_dir = os.path.join(strain_path, f'Layer_{i+1}')
        os.makedirs(layer_dir, exist_ok=True)
        
        # Copy and rename position file
        shutil.copy(pos_files[i], os.path.join(layer_dir, 'pos'))
        print(f"  - Created subdirectory {layer_dir} and copied position file.")
    
    print("\nDirectory setup complete.")


if __name__ == '__main__':
    setup_analysis_dirs()
