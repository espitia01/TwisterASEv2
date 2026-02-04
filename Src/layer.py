from ase import Atoms
import numpy as np
from ase.build import make_supercell


def ensure_proper_transformation(matrix):
    """Ensure transformation matrix is proper (positive determinant)."""
    matrix = np.array(matrix, dtype=float)
    det = np.linalg.det(matrix)
    if abs(det) < 1e-10:
        raise ValueError("Transformation matrix is singular")
    if det < 0:
        matrix[0] *= -1
    matrix = np.round(matrix)
    return matrix


class Layer:
    """
    Represents a single layer in a twisted structure.
    Handles both hexagonal and orthorhombic unit cells.
    """
    
    def __init__(self, layer_file):
        from file_io import parse_layer_file
        
        self.layer_data = parse_layer_file(layer_file)
        self.unitcell = None
        self.unitcell_rot = None
        self.supercell = None
        self.unitcell_ortho = None
        self.unitcell_ortho_rot = None
        self.supercell_ortho = None
        self.has_orthocell = False
        
        self._build_unitcell()
        self._rotate_unitcell()
        
    def _build_unitcell(self):
        """Build the hexagonal unit cell from layer data."""
        data = self.layer_data
        
        lattice_vectors = data['lattice_vectors'] * data['lattice_parameters']
        
        self.unitcell = Atoms(
            symbols=data['atom_symbols'],
            scaled_positions=data['atom_positions'],
            cell=lattice_vectors,
            pbc=[True, True, True]
        )
        
        if data['atom_tags'] is not None:
            self.unitcell.set_tags(data['atom_tags'])
        
        self.unitcell.translate([0, 0, data['translate_z']])
        
        if data['atom_positions_ortho'] is not None:
            self.has_orthocell = True
            orthocell_tm = data['orthocell_tm']
            if orthocell_tm is None:
                orthocell_tm = np.array([[2, 0, 0], [-1, 2, 0], [0, 0, 1]])
            
            ortho_vectors = np.dot(orthocell_tm, lattice_vectors)
            
            self.unitcell_ortho = Atoms(
                symbols=data['atom_symbols_ortho'],
                scaled_positions=data['atom_positions_ortho'],
                cell=ortho_vectors,
                pbc=[True, True, True]
            )
            
            if data['atom_tags_ortho'] is not None:
                self.unitcell_ortho.set_tags(data['atom_tags_ortho'])
            
            self.unitcell_ortho.translate([0, 0, data['translate_z']])
    
    def _rotate_unitcell(self):
        """Apply twist angle rotation to unit cell."""
        twist_angle = self.layer_data['twist_angle']
        
        self.unitcell_rot = self.unitcell.copy()
        self.unitcell_rot.rotate(twist_angle, 'z', rotate_cell=True, center=(0, 0, 0))
        
        if self.has_orthocell:
            self.unitcell_ortho_rot = self.unitcell_ortho.copy()
            self.unitcell_ortho_rot.rotate(twist_angle, 'z', rotate_cell=True, center=(0, 0, 0))
    
    def set_twist_angle(self, new_angle):
        """Override the twist angle and re-rotate the unit cell."""
        self.layer_data['twist_angle'] = new_angle
        self._rotate_unitcell()
    
    def create_supercell(self, sc_matrix_cart):
        """Create supercell from Cartesian supercell matrix."""
        tmatrix = np.dot(sc_matrix_cart, np.linalg.inv(self.unitcell_rot.get_cell()))
        tmatrix = ensure_proper_transformation(tmatrix)
        self.supercell = make_supercell(self.unitcell_rot, tmatrix)
    
    def create_supercell_ortho(self, sc_matrix_cart):
        """Create orthorhombic supercell from Cartesian supercell matrix."""
        if not self.has_orthocell:
            return
        tmatrix = np.dot(sc_matrix_cart, np.linalg.inv(self.unitcell_ortho_rot.get_cell()))
        tmatrix = ensure_proper_transformation(tmatrix)
        self.supercell_ortho = make_supercell(self.unitcell_ortho_rot, tmatrix)
    
    def get_intralayer_potential(self):
        """Return the intralayer potential name."""
        return self.layer_data['intralayer_potential']
