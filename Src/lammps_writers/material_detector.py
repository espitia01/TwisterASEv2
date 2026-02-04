import numpy as np

def detect_material_type(atoms_or_symbols):
    """
    Detect material type from ASE Atoms object or list of symbols.
    
    Returns:
        'TMD', 'Graphene', 'hBN', or 'Unknown'
    """
    if hasattr(atoms_or_symbols, 'get_chemical_symbols'):
        symbols = set(atoms_or_symbols.get_chemical_symbols())
    else:
        symbols = set(atoms_or_symbols)
    
    # TMD detection: Contains transition metal + chalcogen
    TMD_METALS = {'Mo', 'W', 'Nb', 'Ta', 'Re', 'Tc'}
    TMD_CHALCOGENS = {'S', 'Se', 'Te'}
    
    has_metal = bool(symbols & TMD_METALS)
    has_chalcogen = bool(symbols & TMD_CHALCOGENS)
    
    if has_metal and has_chalcogen:
        return 'TMD'
    
    # Graphene detection: Only carbon
    if symbols == {'C'}:
        return 'Graphene'
    
    # hBN detection: Boron and Nitrogen
    if symbols == {'B', 'N'}:
        return 'hBN'
    
    return 'Unknown'


def detect_system_materials(layers):
    """
    Detect all material types present in a list of Layer objects.
    
    Args:
        layers: List of Layer objects with .supercell or .supercell_ortho attributes
    
    Returns:
        dict: {
            'layer_types': list of material types per layer,
            'has_hbn': bool,
            'has_graphene': bool,
            'has_tmd': bool,
            'requires_full_format': bool (True if hBN present)
        }
    """
    layer_types = []
    
    for layer in layers:
        # Use ortho supercell if available, otherwise hex supercell
        supercell = layer.supercell_ortho if layer.has_orthocell and layer.supercell_ortho is not None else layer.supercell
        material_type = detect_material_type(supercell)
        layer_types.append(material_type)
    
    has_hbn = 'hBN' in layer_types
    has_graphene = 'Graphene' in layer_types
    has_tmd = 'TMD' in layer_types
    
    # Use full format (with molecule ID and charge) if hBN is present
    requires_full_format = has_hbn
    
    return {
        'layer_types': layer_types,
        'has_hbn': has_hbn,
        'has_graphene': has_graphene,
        'has_tmd': has_tmd,
        'requires_full_format': requires_full_format
    }
