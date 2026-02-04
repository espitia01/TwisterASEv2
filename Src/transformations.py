import numpy as np
from ase import Atoms
from ase.build import cut
from scipy.spatial import cKDTree


def angle_tm_using_i(i_val):
    """
    Calculate twist angle and transformation matrix using i parameter.
    Reference: https://doi.org/10.1016/j.cpc.2021.108184
    """
    transformation_matrix = [
        [i_val, i_val + 1, 0],
        [-(i_val + 1), 2 * i_val + 1, 0],
        [0, 0, 1]
    ]
    theta = np.arccos((3 * i_val**2 + 3 * i_val + 0.5) / 
                      (3 * i_val**2 + 3 * i_val + 1)) * 180 / np.pi
    return theta, np.array(transformation_matrix)


def angle_tm_using_mn(m_val, n_val):
    """
    Calculate twist angle and transformation matrix using m,n parameters.
    """
    transformation_matrix = [
        [m_val, n_val, 0],
        [m_val + n_val, -1 * m_val, 0],
        [0, 0, 1]
    ]
    theta = np.arccos((m_val**2 + 4 * m_val * n_val + n_val**2) / 
                      (2 * (m_val**2 + m_val * n_val + n_val**2))) * 180. / np.pi
    return theta, np.round(transformation_matrix)


def compute_ortho_supercell_matrix(hex_supercell_vectors, ortho_unitcell_vectors):
    """
    Compute transformation matrix from hexagonal to orthorhombic supercell.
    """
    hex_to_orth = np.copy(hex_supercell_vectors)
    hex_to_orth[0] = hex_supercell_vectors[0]
    hex_to_orth[1] = 2 * hex_supercell_vectors[1] - hex_supercell_vectors[0]
    
    transformation_matrix = np.dot(hex_to_orth, np.linalg.inv(ortho_unitcell_vectors))
    transformation_matrix[0:2] = transformation_matrix[0:2] * 4
    transformation_matrix = transformation_matrix.round().astype(int)
    
    return np.dot(transformation_matrix, ortho_unitcell_vectors)


def remove_overlapping_atoms_ckdtree(positions, tolerance):
    """
    Remove overlapping atoms using cKDTree for efficient neighbor search.
    
    Parameters:
        positions: Array of shape (N, 3) with atom coordinates
        tolerance: Distance below which atoms are considered overlapping
    
    Returns:
        filtered_positions: Array of unique positions
        ind_uniq: Indices of unique atoms
    """
    tree = cKDTree(positions)
    visited = np.zeros(len(positions), dtype=bool)
    filtered_positions = []
    ind_uniq = []
    
    for i, pos in enumerate(positions):
        if not visited[i]:
            filtered_positions.append(pos)
            ind_uniq.append(i)
            pos = pos.copy(order='C')
            indices = tree.query_ball_point(pos, tolerance)
            visited[indices] = True
    
    return np.array(filtered_positions), np.array(ind_uniq)


def fold(rcrys):
    """Fold crystal coordinates into range [0, 1)."""
    rfold = np.copy(rcrys)
    for iat in range(len(rcrys)):
        for jj in range(3):
            while rfold[iat, jj] < -1e-5:
                rfold[iat, jj] = rfold[iat, jj] + 1.0
            while rfold[iat, jj] >= 1 - 1e-5:
                rfold[iat, jj] = rfold[iat, jj] - 1.0
    return rfold


def hexcut(struct, cut_to_lattice, cut_from_lattice):
    """
    Cut structure to a smaller supercell.
    
    Parameters:
        struct: ASE Atoms object
        cut_to_lattice: Target lattice vectors (3x3 array)
        cut_from_lattice: Source lattice vectors (3x3 array)
    
    Returns:
        struct_cutsl: ASE Atoms object with cut structure
    """
    tr = np.dot(cut_to_lattice, np.linalg.inv(cut_from_lattice))
    cut_sl = cut(struct, a=tr[0], b=tr[1], c=tr[2], tolerance=0.3, extend=1.2)
    
    cut_sl_pos_c = cut_sl.get_scaled_positions()
    cut_sl_sym = np.array(cut_sl.get_chemical_symbols())
    cut_sl_tags = cut_sl.get_tags()
    cut_sl_pos_c = fold(cut_sl_pos_c)
    
    cut_sl_pos_c, indices = remove_overlapping_atoms_ckdtree(cut_sl_pos_c, 1e-3)
    
    struct_cutsl = Atoms(
        cell=cut_to_lattice,
        symbols=cut_sl_sym[indices],
        scaled_positions=cut_sl_pos_c,
        pbc=[True, True, True]
    )
    struct_cutsl.set_tags(cut_sl_tags[indices])
    
    return struct_cutsl
