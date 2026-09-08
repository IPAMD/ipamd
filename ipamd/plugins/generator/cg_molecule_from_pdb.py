from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.parser import range_to_list
import os

configure = {
    "resource": ['persistency_dir', 'ff'],
}
def func(file_name, persistency_dir=None, ff=None, rigid_range=''):
    position_of_last_dot = file_name.rfind('.')
    protein_name = file_name[:position_of_last_dot]
    molecule = Molecule(protein_name, cg='CM')
    rigid_atoms_index = []
    if rigid_range != '':
        tmp = range_to_list(rigid_range)
        for i in tmp:
            rigid_atoms_index.append(i - 1)
    rigid_groups = []
    for i in rigid_atoms_index:
        if len(rigid_groups) == 0:
            rigid_groups.append([i])
        else:
            if i - rigid_groups[-1][-1] <= 1:
                rigid_groups[-1].append(i)
            else:
                rigid_groups.append([i])
    rigid_groups_filtered = []
    for group in rigid_groups:
        if len(group) >= 2:
            rigid_groups_filtered.append(group)
    rigid_groups = rigid_groups_filtered
    full_path = os.path.join(persistency_dir, file_name)
    with open(full_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        residue_id = 0
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                charge = line[78:80].strip()

                rigid_group = -1
                for group_id, group in enumerate(rigid_groups):
                    if residue_id in group:
                        rigid_group = group_id
                        break
                molecule.add_atom(
                    atom=Atom(
                        velocity=[0, 0, 0],
                        atom_type=name,
                        ff = ff,
                        charge=float(charge) if charge != '' else None
                    ),
                    coordinate=[x / 10, y / 10, z / 10],
                    rigid_group=rigid_group,
                )
            elif line.startswith('CONNECT'):
                items = line.split()
                atom1 = int(items[1]) - 1
                atom2 = int(items[2]) - 1
                molecule.link(atom1, atom2, 'B-B')
    return molecule
