from ipamd.public.models.md import Molecule, Atom
from ipamd.public.utils.math import polar_to_cartesian
from ipamd.public.models.sequence import ProteinSequence
import math
configure = {
    "resource": ['ff'],
}
def func(protein: ProteinSequence, ff=None, d=0.38):
    molecule = Molecule(protein.name, cg='CA')
    length = len(protein)
    r = 0.38
    b = d / (2 * math.pi)
    phi = r / b
    for index, res in enumerate(protein.sequence):
        if index == 0:
            custom_mass = f'compute:{ff.atom_definition[res]["mass"]}+1.008'
        elif index == length - 1:
            custom_mass = f'compute:{ff.atom_definition[res]["mass"]}+17.007'
        else:
            custom_mass = None
        x, y = polar_to_cartesian(r, phi)
        molecule.add_atom(
            Atom(
                velocity=(0, 0, 0),
                atom_type=res,
                mass=custom_mass,
                ff=ff
            ),
            coordinate=(x, y, 0),
            bond='B-B'
        )
        phi += 0.38 / r
        r = b * phi
    return molecule
