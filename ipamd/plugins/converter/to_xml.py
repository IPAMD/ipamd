import os
import math
from ipamd.public.utils.xml import dict_to_xml
configure = {
    'type': 'function',
    "schema": 'frame',
    "apply": ['persistency_dir']
}
def func(filename, persistency_dir, frame, ignoring_pbc = False, indent = True):
    prop = frame.properties(ignoring_image=ignoring_pbc)
    conf_path = os.path.join(persistency_dir, filename + '.xml')
    map_path = os.path.join(persistency_dir, filename + '.map')
    extras = []
    content_type = ''
    content_position = ''
    content_velocity = ''
    content_image = ''
    content_mass = ''
    content_charge = ''
    content_bond = ''
    content_molecule = ''
    content_body = ''
    n_bonds = 0
    n_atoms = 0
    n_molecule = 0

    epsilon = frame.box.env.values['epsilon']
    total_rigid_groups = 0
    for molecule in prop['molecules']:
        extras.append(
            {
                'type': molecule['molecule_type'],
                'CG': molecule['CG']
            }
        )
        content_type += '\n' .join(molecule['type']) + '\n'
        content_position += '\n' .join(['\t'.join(map(str, pos)) for pos in molecule['position']]) + '\n'
        content_velocity += '\n'.join(['\t'.join(map(str, v)) for v in molecule['velocity']]) + '\n'
        content_image += '\n' .join(['\t'.join(map(lambda s: str(int(s)), img)) for img in molecule['image']]) + '\n'
        content_mass += '\n' .join(map(str, molecule['mass'])) + '\n'
        content_charge += '\n' .join(map(lambda x: str(x * math.sqrt(138.935 / epsilon)), molecule['charge'])) + '\n'
        for i in range(molecule['n_atoms']):
            content_molecule += str(n_molecule) + '\n'
        content_body += '\n' .join(map(lambda x: str(-1 if x == -1 else x + total_rigid_groups), molecule['rigid_group'])) + '\n'
        max_group_index = max(molecule['rigid_group'])
        total_rigid_groups += max_group_index + 1
        content_bond += '\n'.join(
            '\t'.join(
                map(str,
                    [bond['type'], bond['atom1'] + n_atoms, bond['atom2'] + n_atoms]
                )
            )
            for bond in molecule['bonds']
        ) + '\n'
        n_atoms += molecule['n_atoms']
        n_bonds += molecule['n_bonds']
        n_molecule += 1
    d = {
        'galamost_xml': {
            'version': 1.3,
            'configuration': {
                'time_step': "0",
                'dimensions': "3",
                'natoms': n_atoms,
                'box': {
                    'lx': prop['size'][0],
                    'ly': prop['size'][1],
                    'lz': prop['size'][2]
                },
                'type': {
                    'num': n_atoms,
                    'text': content_type
                },
                'body': {
                    'num': n_atoms,
                    'text': content_body
                },
                'position': {
                    'num': n_atoms,
                    'text': content_position
                },
                'velocity': {
                    'num': n_atoms,
                    'text': content_velocity
                },
                'image': {
                    'num': n_atoms,
                    'text': content_image
                },
                'mass': {
                    'num': n_atoms,
                    'text': content_mass
                },
                'charge': {
                    'num': n_atoms,
                    'text': content_charge
                },
                'molecule': {
                    'num': n_atoms,
                    'text': content_molecule
                },
                'bond': {
                    'num': n_bonds,
                    'text': content_bond
                }

            }
        }
    }
    with open(conf_path, 'w', encoding='utf-8') as f:
        f.write(dict_to_xml(d, indent=indent))
    with open(map_path, 'w', encoding='utf-8') as f:
        for i, extra in enumerate(extras):
            f.write(f'{i} {extra['type']} {extra['CG']}\n')
