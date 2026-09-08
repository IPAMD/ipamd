from ipamd.public.models.md import Molecule, Atom
from ipamd.public.models.sequence import ProteinSequence
import math
import copy
configure = {
    "resource": ['ff'],
}
def func(protein: ProteinSequence, ff=None, span=2):
    molecule = Molecule(protein.name, cg='CA')
    length = len(protein)
    n_group = math.ceil(length / span)
    width = math.ceil(math.pow(n_group, 1 / 3))
    position = [0.0, 0.0, 0.0]
    for group_index in range(n_group):
        plane_index, index_in_plane = divmod(group_index, width * width)
        row_index, index_in_row = divmod(group_index, width)

        if index_in_plane + 1 == width * width:
            dimension = 2
            direction = 1
        elif index_in_row + 1 == width:
            dimension = 1
            direction = 1 if (plane_index % 2 == 0) else -1
        else:
            dimension = 0
            direction = 1 if (row_index % 2 == 0) else -1

        for i in range(
            min(span, length - group_index * span)
        ):
            index = group_index * span + i
            res = protein.sequence[index]
            if index == 0:
                custom_mass = f'compute:{ff.atom_definition[res]["mass"]}+1.008'
            elif index == length - 1:
                custom_mass = f'compute:{ff.atom_definition[res]["mass"]}+17.007'
            else:
                custom_mass = None

            molecule.add_atom(
                Atom(
                    velocity=(0, 0, 0),
                    atom_type=res,
                    mass=custom_mass,
                    ff=ff
                ),
                coordinate=copy.deepcopy(position),
                bond='B-B'
            )
            position[dimension] += 0.38 * direction
    return molecule