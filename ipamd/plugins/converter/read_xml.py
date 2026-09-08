import os.path
from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.output import warning
from ipamd.public.utils.parser import value_of
from ipamd.public.utils.xml import read_xml
import re
import math

configure = {
    'type': 'function',
    "schema": 'frame',
    "apply": ['persistency_dir', 'ff']
}
def func(filename, persistency_dir, frame, ff):
    xml_path = os.path.join(persistency_dir, filename + '.xml')
    if len(re.findall(r'\.\d+', filename)) == 0:
        filename = filename
    else:
        position_of_dot = filename.rfind('.')
        filename = filename[:position_of_dot]
    map_path = os.path.join(persistency_dir, filename + '.map')
    molecule_map = []
    cg_map = []
    env = frame.box.env.values
    try:
        with open(map_path, 'r', encoding='utf-8') as f:
            map_content = f.read()
            map_content = map_content.split('\n')
            for line in map_content:
                if line == '':
                    continue
                else:
                    _, molecule_name, cg = line.split()
                    molecule_map.append(molecule_name)
                    cg_map.append(cg)
    except FileNotFoundError:
        warning(f'molecule name map({map_path}) not found, all molecules will be named mol')

    content = read_xml(xml_path)
    contained_molecules = []
    molecule_configure = content['galamost_xml']['configuration']
    size = molecule_configure['box']
    x = float(size['lx'])
    y = float(size['ly'])
    z = float(size['lz'])
    frame.set_size(x, y, z)
    atom_num = int(molecule_configure['type']['num'])

    def get_content(name):
        if name in molecule_configure.keys():
            node_txt = molecule_configure[name]['text']
        else:
            node_txt = ''
        return list(
            map(
                lambda s: s.strip(),
                node_txt.split('\n')
            )
        )

    type_content = get_content('type')
    body_content =  get_content('body')
    mass_content = get_content('mass')
    charge_content = get_content('charge')
    position_content = get_content('position')
    image_content = get_content('image')
    molecule_index_content = get_content('molecule')
    bond_content = get_content('bond')
    velocity_content = get_content('velocity')
    angle_content = get_content('angle')

    current_molecule_index = 0
    current_atom_index_in_mol = 0
    mn = molecule_map[current_molecule_index] if len(molecule_map) != 0 else 'mol'
    cg = cg_map[current_molecule_index] if len(cg_map) != 0 else 'CA'
    molecule = Molecule(mn, cg)
    atom_map = []
    previous_rigid_group = 0
    for current_atom_index in range(atom_num):
        if current_molecule_index != int(molecule_index_content[current_atom_index]):
            current_molecule_index += 1
            body_prop = molecule.properties()['rigid_group']
            previous_rigid_group += max(body_prop) + 1
            contained_molecules.append(molecule)
            mn = molecule_map[current_molecule_index] if len(molecule_map) != 0 else 'mol'
            cg = cg_map[current_molecule_index] if len(cg_map) != 0 else 'CA'
            molecule = Molecule(mn, cg)
            current_atom_index_in_mol = 0


        atom_map.append((current_molecule_index, current_atom_index_in_mol))
        current_atom_index_in_mol += 1
        rigid_group = int(body_content[current_atom_index]) - previous_rigid_group
        position = list(map(float, position_content[current_atom_index].split()))
        velocity = list(map(float, velocity_content[current_atom_index].split()))
        image = list(map(int, image_content[current_atom_index].split()))
        charge = float(charge_content[current_atom_index]) / math.sqrt(138.935 / frame.box.env.values['epsilon'])
        mass = float(mass_content[current_atom_index])
        atom_type = type_content[current_atom_index]

        mass_excepted = value_of(ff.atom_definition[atom_type]['mass'], env)
        charge_excepted = value_of(ff.atom_definition[atom_type]['charge'], env)
        if abs(mass_excepted - mass) < 1e-6:
            mass = None
        if abs(charge_excepted - charge) < 1e-6:
            charge = None

        position = [
            position[0] + image[0] * x + 0.5 * x,
            position[1] + image[1] * y + 0.5 * y,
            position[2] + image[2] * z + 0.5 * z
        ]
        atom = Atom(
            velocity=velocity,
            atom_type=atom_type,
            mass=mass,
            charge=charge,
            ff=ff
        )
        molecule.add_atom(
            atom,
            coordinate=position,
            rigid_group=rigid_group,
            bond=None
        )
    contained_molecules.append(molecule)

    for bond in bond_content:
        if bond == '':
            continue
        bond_type, index1, index2 = bond.split()
        index1 = int(index1)
        index2 = int(index2)
        molecule_id = atom_map[index1][0]
        atom_id1 = atom_map[index1][1]
        atom_id2 = atom_map[index2][1]
        contained_molecules[molecule_id].link(atom_id1, atom_id2, type_=bond_type)

    for angle in angle_content:
        if angle == '':
            continue
        angle_type, index1, index2, index3 = angle.split()
        index1 = int(index1)
        index2 = int(index2)
        index3 = int(index3)
        molecule_id = atom_map[index1][0]
        atom_id1 = atom_map[index1][1]
        atom_id2 = atom_map[index2][1]
        atom_id3 = atom_map[index3][1]
        contained_molecules[molecule_id].add_angle(atom_id1, atom_id2, atom_id3, type_=angle_type)

    for molecule in contained_molecules:
        frame.add_molecule(molecule)
