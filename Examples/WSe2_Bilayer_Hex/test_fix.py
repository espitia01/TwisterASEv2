#!/usr/bin/env python3
"""Test the TMD input generator fix"""

# Simulate the fixed code
material_info = {
    'atom_types': [1, 2, 3, 4, 5, 6],
    'type_to_symbol': {
        1: 'W',
        2: 'Se', 
        3: 'Se',
        4: 'W',
        5: 'Se',
        6: 'Se'
    },
    'layer_info': [
        {'unique_tags': [1, 2, 3], 'intralayer_potential': 'tmd.sw'},
        {'unique_tags': [4, 5, 6], 'intralayer_potential': 'tmd.sw'}
    ]
}

print("Testing fixed TMD intralayer interaction generation:\n")

for i, layer_info in enumerate(material_info['layer_info']):
    # Create symbol mapping for all atom types
    symbols_line = ['NULL'] * len(material_info['atom_types'])
    # Map this layer's tags to their symbols
    for tag in layer_info['unique_tags']:
        symbols_line[tag - 1] = material_info['type_to_symbol'][tag]
    pot = layer_info['intralayer_potential']
    print(f"Layer {i+1}:")
    print(f"pair_coeff * * sw/mod {i+1} {pot} {' '.join(symbols_line)}")
    print()

print("\nExpected output:")
print("Layer 1: W Se Se NULL NULL NULL")
print("Layer 2: NULL NULL NULL W Se Se")
