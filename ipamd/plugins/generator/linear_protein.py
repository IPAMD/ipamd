from ipamd.public.models.md import Molecule, Atom
from ipamd.public.models.sequence import ProteinSequence

configure = {
    "resource": ['ff'],
}
def func(protein: ProteinSequence, ff=None):
    molecule = Molecule(protein.name, cg='CA')
    length = len(protein)
    z = 0.0
    for index, r in enumerate(protein.sequence):
        if index == 0:
            custom_mass = f'compute:{ff.atom_definition[r]["mass"]}+1.008'
        elif index == length - 1:
            custom_mass = f'compute:{ff.atom_definition[r]["mass"]}+17.007'
        else:
            custom_mass = None
        molecule.add_atom(
            Atom(
                velocity=(0, 0, 0),
                atom_type=r,
                mass=custom_mass,
                ff=ff
            ),
            coordinate=(0.0, 0.0, z),
            bond='B-B'
        )
        z += 0.38
    return molecule
