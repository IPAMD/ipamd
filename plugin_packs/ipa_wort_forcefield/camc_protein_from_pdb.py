import os
import periodictable
from pypdbio import PdbReader
from pybioseq.converter import standard_aa, to1
from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.output import warning, info, error

configure = {
    "resource": ['persistency_dir', 'ff'],
}

def func(file_name, persistency_dir=None, ff=None):
    reader = PdbReader(os.path.join(persistency_dir, file_name))
    res = reader.read()
    if res.meta.pdb_id:
        protein_name = res.meta.pdb_id
    else:
        protein_name, _ = os.path.splitext(file_name)

    molecule = Molecule(protein_name, cg='CAMC')
    model = res.models[0]

    n_chains = len(model.chains)
    if n_chains > 1:
        warning('more than one chain detected.')

    sequence = ''
    for chain in model:
        last_ca_index = -1
        current_index = -1
        for residue in chain:
            residue_name = residue.name
            if standard_aa(residue_name):
                sequence += to1(residue_name)
            else:
                if residue_name == 'HOH':
                    continue
                error('meet a unknown residue: ' + residue_name)
                raise ValueError

            total_mass = 0
            total_mass_x = 0
            total_mass_y = 0
            total_mass_z = 0
            for atom in residue:
                if atom.name == 'CA':
                    ca_coordinate = [atom.coord[0], atom.coord[1], atom.coord[2]]
                    current_index += 1
                    molecule.add_atom(
                        atom=Atom(
                            velocity=[0, 0, 0],
                            atom_type="BB",
                            ff=ff
                        ),
                        coordinate=ca_coordinate,
                    )
                    if last_ca_index != -1:
                        molecule.link(last_ca_index, current_index, 'BB-BB')
                    last_ca_index = current_index
                    continue
                if atom.name in ('C', 'N', 'O'):
                    continue
                if atom.element == 'H':
                    continue
                atom_mass = periodictable.elements.symbol(atom.element).mass
                total_mass += atom_mass
                total_mass_x += atom_mass * atom.coord[0]
                total_mass_y += atom_mass * atom.coord[1]
                total_mass_z += atom_mass * atom.coord[2]

            if total_mass > 0:
                cs_coordinate = [
                    total_mass_x / total_mass,
                    total_mass_y / total_mass,
                    total_mass_z / total_mass
                ]
                molecule.add_atom(
                    atom=Atom(
                        velocity=[0, 0, 0],
                        atom_type=to1(residue_name) + 'S',
                        ff=ff
                    ),
                    coordinate=cs_coordinate,
                )
                current_index += 1
                molecule.link(last_ca_index, current_index, 'B-S')

    info('processed protein with sequence: ' + sequence)
    return molecule
