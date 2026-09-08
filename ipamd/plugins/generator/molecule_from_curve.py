from ipamd.public.models.md import Molecule, Atom
from ipamd.public.models.sequence import Sequence

configure = {
    "resource": ['ff'],
}

def func(
    mol: Sequence, curve, start_add=0, end_add=0, cg="MC", rename_map={}, new_bond="B-B", ff=None
):
    molecule = Molecule(mol.name, cg=cg)
    for index, r in enumerate(mol.sequence):
        remapped_type = rename_map.get(r, r)
        if index == 0:
            custom_mass = f'compute:{ff.atom_definition[remapped_type]["mass"]}+{start_add}'
        elif index == len(mol) - 1:
            custom_mass = f'compute:{ff.atom_definition[remapped_type]["mass"]}+{end_add}'
        else:
            custom_mass = None

        x, y, z = curve(index)
        molecule.add_atom(
            Atom(
                velocity=(0, 0, 0),
                atom_type=remapped_type,
                mass=custom_mass,
                ff=ff
            ),
            coordinate=(x, y, z),
            bond=new_bond
        )
    return molecule
