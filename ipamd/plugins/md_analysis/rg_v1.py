"""
plugin for calculating radius of gyration
"""
import math
from ipamd.public.models.data import Scalar
from ipamd.public.utils.parser import protein_range_split

configure = {
    "schema": ['frame'],
}

def func(frame, target_molecule=None, **kwargs):
    """
    main function for calculating radius of gyration

    :param target_molecule: the target molecule to calculate. If None, calculate for all molecules.
    :param kwargs: other parameters, don't need to be set

    :return: average radius of gyration over selected molecules
    """
    prop = frame.properties(ignoring_image=True)
    molecules = prop['molecules']
    n_molecules = 0
    sum_rg = 0

    target_molecule_type, target_range = protein_range_split(target_molecule)
    for molecule in molecules:
        molecule_type = molecule['molecule_type']
        if target_molecule_type != "" and molecule_type != target_molecule_type:
            continue
        position = molecule['position']
        mass = molecule['mass']
        molecule_length = len(mass)

        sum_mass = 0
        com = [0.0, 0.0, 0.0]
        for i in range(molecule_length):
            if target_range != [] and i not in target_range:
                continue
            sum_mass += mass[i]
            com[0] += mass[i] * position[i][0]
            com[1] += mass[i] * position[i][1]
            com[2] += mass[i] * position[i][2]
        if sum_mass == 0:
            continue
        com = [c / sum_mass for c in com]

        rg2 = 0.0
        for i in range(molecule_length):
            if target_range != [] and i not in target_range:
                continue
            dx = position[i][0] - com[0]
            dy = position[i][1] - com[1]
            dz = position[i][2] - com[2]
            rg2 += mass[i] * (dx * dx + dy * dy + dz * dz)
        sum_rg += math.sqrt(rg2 / sum_mass)
        n_molecules += 1

    return Scalar(
        data=sum_rg / n_molecules if n_molecules > 0 else 0,
        title='rg',
        unit='nm'
    )
