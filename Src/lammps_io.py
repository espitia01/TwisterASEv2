import numpy as np
import os
from ase import Atoms
from ase.io import write

def write_lammps_data_file(filename, ase_atom_struct):
    old_lattice_vectors = np.array(ase_atom_struct.get_cell())
    lattice_vectors = restricted_triclinic(old_lattice_vectors)
    pos_new = write2lammps(filename, ase_atom_struct, lattice_vectors)
    struct_new = Atoms(
    	symbols=ase_atom_struct.symbols,
    	positions=pos_new,
   		cell=lattice_vectors,
    	pbc=[True, True, True]  )
    write('superlattice_lammps.cif', struct_new)
    print(f"Written superlattice_lammps.cif with {len(struct_new)} atoms")

def restricted_triclinic(lattice_vectors):
    '''
    Converts lattice_vectors to a restricted triclinic cell
    '''

    a, b, c = np.array(lattice_vectors[0]), np.array(lattice_vectors[1]), np.array(lattice_vectors[2])
    
    # Step 1: Compute ê₁ (normalized a)
    e1 = a / np.linalg.norm(a)
    
    # Step 2: Compute ê₂ (orthogonalized b)
    b_proj = np.dot(b, e1) * e1
    b_perp = b - b_proj
    if np.linalg.norm(b_perp) < 1e-10:
        raise ValueError("a and b are colinear; no unique basis")
    e2 = b_perp / np.linalg.norm(b_perp)
    
    # Step 3: Compute ê₃ (cross product)
    e3 = np.cross(e1, e2)
    
    # Step 4: Verify orientation constraint for c'
    if np.dot(e3, c) <= 0:
        raise ValueError("c' would have non-positive z-component")
    
    # Construct change-of-basis matrix
    O = np.vstack([e1, e2, e3])
    return (O @ lattice_vectors.T).T

def write2lammps(filename, atoms, lattice_vectors):
    with open (filename, 'w') as f:
        f.write('LAMMPS data file\n\n')
        # Number of atoms
        f.write('{} atoms\n'.format(len(atoms.positions)))
        # Number of atom types are number of  unique elements in atoms.tags
        f.write('{} atom types\n'.format(len(np.unique(atoms.get_tags()))))
        # Box bounds from lattice_vectors
        f.write('0.0 {} xlo xhi\n'.format(lattice_vectors[0,0]))
        f.write('0.0 {} ylo yhi\n'.format(lattice_vectors[1,1]))
        f.write('0.0 {} zlo zhi\n'.format(lattice_vectors[2,2]))
        f.write('{} {} {} xy xz yz\n'.format(lattice_vectors[1,0], lattice_vectors[2,0], lattice_vectors[2,1]))
        f.write('\n')
        f.write('Masses\n\n')
        # find unique tags in atoms.get_tags() and their indices
        unique_tags, indices = np.unique(atoms.get_tags(), return_index=True)
        for ind in indices:
            f.write('{} {}\n'.format(atoms.get_tags()[ind], atoms.get_masses()[ind])) 
        f.write('\n')

        f.write('Atoms\n\n')
        tags = atoms.get_tags()
        pos_arr = np.zeros((len(atoms.positions), 3))
        for i, iat in enumerate(atoms.get_scaled_positions()):
            # convert scaled positions to Cartesian
            pos = np.dot(iat, lattice_vectors)
            pos_arr[i] = pos
            f.write('{} {} {} {} {}\n'.format(i+1, tags[i], pos[0], pos[1], pos[2]))
    return pos_arr
