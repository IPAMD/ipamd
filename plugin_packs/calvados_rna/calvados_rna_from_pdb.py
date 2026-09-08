import os
from pypdbio import PdbReader
from pybioseq.converter import standard_nmp
from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.output import warning, info

configure = {
    "resource": ['persistency_dir', 'ff'],
}

def func(file_name, persistency_dir=None, ff=None):
    reader = PdbReader(os.path.join(persistency_dir, file_name))
    res = reader.read()
    if res.meta.pdb_id:
        rna_name = res.meta.pdb_id
    else:
        rna_name, _ = os.path.splitext(file_name)

    molecule = Molecule(rna_name, cg='CARNA')
    model = res.models[0]

    if len(model.chains) > 1:
        warning('more than one chain detected.')

    sequence = ''
    for chain in model:
        for index, residue in enumerate(chain):
            residue_name = residue.name
            if standard_nmp(residue_name):
                sequence += residue_name
            else:
                if residue_name == 'HOH':
                    continue
                warning('meet a unknown residue: ' + residue_name)
                continue
            atom_p = residue["P"]
            atom_n = None
            if residue_name == "U" or residue_name == "A" or residue_name == "C" :
                atom_n = residue["N1"]
            elif residue_name == "G":
                atom_n = residue["N9"]

            molecule.add_atom(
                atom=Atom(
                    velocity=[0, 0, 0],
                    atom_type='NMPB',
                    ff=ff
                ),
                coordinate=[atom_p.coord[0], atom_p.coord[1], atom_p.coord[2]],
                bond=None
            )
            molecule.add_atom(
                atom=Atom(
                    velocity=[0, 0, 0],
                    atom_type=f'NMP{residue_name}',
                    ff=ff
                ),
                coordinate=[atom_n.coord[0], atom_n.coord[1], atom_n.coord[2]],
                bond=None
            )
            molecule.link(2 * index, 2 * index + 1, 'RB-RS')
            if index > 0:
                molecule.link(
                    2 * index,
                    2 * index - 2,
                    'RB-RB'
                )
                molecule.link(2 * index + 1, 2 * index - 1, 'RS-RS')

        chain_length = len(sequence)
        for i in range(1, chain_length - 1):
            molecule.add_angle(2 * i - 2, 2 * i, 2 * i + 2, 'RB-RB-RB')

    info('processed RNA with sequence: ' + sequence)
    return molecule
