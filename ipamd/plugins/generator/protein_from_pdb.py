import os
import periodictable
from pypdbio import PdbReader
from pybioseq.converter import standard_aa, to1
from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.output import warning, info, error
from ipamd.public.utils.parser import range_to_list

configure = {
    "resource": ['persistency_dir', 'ff'],
}
def func(
    file_name,
    persistency_dir=None,
    ff=None,
    ignoring_h=True,
    cg='mass_center',
    rigid_range='',
    rigid_from_plddt=False,
    threshold=70,
    max_gap=4
):
    reader = PdbReader(os.path.join(persistency_dir, file_name))
    res = reader.read()
    if res.meta.pdb_id:
        protein_name = res.meta.pdb_id
    else:
        protein_name, _ = os.path.splitext(file_name)
    cg = 'CA' if cg == 'alpha' else 'MC'
    molecule = Molecule(protein_name, cg=cg)
    model = res.models[0]
    #determine the rigid atoms
    rigid_groups = []
    if rigid_range != '':
        groups = []
        if ';' in rigid_range:
            groups = rigid_range.split(';')
        else:
            groups = [rigid_range]
        for group_str in groups:
            rigid_groups.append(range_to_list(group_str))
    elif rigid_from_plddt:
        for chain in model:
            for residue_id, residue in enumerate(chain.residues):
                plddt = residue.atoms[0].temp_factor * 100
                if plddt > threshold:
                    if len(rigid_groups) == 0:
                        rigid_groups.append([residue_id])
                    else:
                        if residue_id - rigid_groups[-1][-1] <= max_gap + 1:
                            rigid_groups[-1].append(residue_id)
                        else:
                            rigid_groups.append([residue_id])

    rigid_groups_filtered = []
    for group in rigid_groups:
        if len(group) >= 2:
            rigid_groups_filtered.append(group)
    rigid_groups = rigid_groups_filtered

    n_chains = len(model.chains)
    if n_chains > 1:
        warning('more than one chain detected.')

    residue_id = 0
    for chain in model:
        new_chain = True
        sequence = ''
        chain_length = len(chain.residues)
        for residue in chain:
            residue_name = residue.name
            if standard_aa(residue_name):
                residue_id += 1
                sequence += to1(residue_name)
            else:
                if residue_name == 'HOH':
                    continue
                else:
                    error('meet a unknown residue: ' + residue_name)
                    raise ValueError
            ca_coordinate = [0, 0, 0]
            mc_coordinate = [0, 0, 0]
            total_mass = 0
            for atom in residue:
                atom_mass = periodictable.elements.symbol(atom.element).mass
                if atom.name == 'CA':
                    ca_coordinate = [atom.coord[0], atom.coord[1], atom.coord[2]]
                if atom.element != 'H':
                    total_mass += atom_mass
                    mc_coordinate[0] += atom_mass * atom.coord[0]
                    mc_coordinate[1] += atom_mass * atom.coord[1]
                    mc_coordinate[2] += atom_mass * atom.coord[2]
                else:
                    if not ignoring_h:
                        total_mass += atom_mass
                        mc_coordinate[0] += atom_mass * atom.coord[0]
                        mc_coordinate[1] += atom_mass * atom.coord[1]
                        mc_coordinate[2] += atom_mass * atom.coord[2]
            mc_coordinate[0] /= total_mass
            mc_coordinate[1] /= total_mass
            mc_coordinate[2] /= total_mass
            rigid_group = -1
            for group_id, group in enumerate(rigid_groups):
                if residue_id in group:
                    rigid_group = group_id
                    break
            if rigid_group != -1 and cg == 'MC':
                coordinate = mc_coordinate
            else:
                coordinate = ca_coordinate
            res_name_one_word = to1(residue_name)
            if residue_id == 0:
                custom_mass = f'compute:{ff.atom_definition[res_name_one_word]["mass"]}+1.008'
            elif residue_id == chain_length - 1:
                custom_mass = f'compute:{ff.atom_definition[res_name_one_word]["mass"]}+17.007'
            else:
                custom_mass = None
            molecule.add_atom(
                atom=Atom(
                    velocity=[0, 0, 0],
                    atom_type=res_name_one_word,
                    mass=custom_mass,
                    ff = ff
                ),
                coordinate=coordinate,
                rigid_group=rigid_group,
                bond=None if new_chain else 'B-B'
            )
            new_chain = False
    info('processed protein with sequence: ' + sequence)
    return molecule
