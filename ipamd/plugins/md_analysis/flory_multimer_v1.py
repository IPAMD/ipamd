"""
plugin for estimating the Flory scaling exponent from chains of mixed lengths
"""
import math
import numpy as np
from ipamd.public.models.data import Scalar
from ipamd.public.utils.output import warning, error
from ipamd.public.utils.parser import frame_index_list, protein_range_split

def func(box, target_frame, target_molecule=None, **kwargs):
    indices = frame_index_list(target_frame)
    if len(indices) == 0:
        indices = [box.current_frame().no]
    if len(indices) <= 5:
        warning("Small number of frames may lead to inaccurate results.")

    target_molecule_type, target_range = protein_range_split(target_molecule)
    n_list = []
    r_list = []
    for index in indices:
        prop = box.frame(index).current_frame().properties(ignoring_image=True)
        for molecule in prop['molecules']:
            if target_molecule_type != "" and molecule['molecule_type'] != target_molecule_type:
                continue
            mass = molecule['mass']
            position = molecule['position']
            if target_range:
                mass = [mass[i] for i in target_range]
                position = [position[i] for i in target_range]
            molecule_length = len(mass)
            total_mass = 0.0
            mass_center = np.zeros(3, dtype=float)
            for i in range(molecule_length):
                total_mass += mass[i]
                mass_center += mass[i] * np.array(position[i], dtype=float)
            mass_center /= total_mass
            r2_m = 0.0
            for i in range(molecule_length):
                distance2 = np.linalg.norm(np.array(position[i], dtype=float) - mass_center) ** 2
                r2_m += distance2 * mass[i]
            rg2 = r2_m / total_mass
            n_list.append(math.log2(molecule_length))
            r_list.append(math.log2(rg2) / 2)

    if len(n_list) == 0:
        error('No molecule matched the selected type')
        return Scalar(data=0, title='Flory exponent', unit='')

    design = np.vstack([n_list, np.ones(len(n_list))]).T
    slope, _ = np.linalg.lstsq(design, r_list, rcond=None)[0]
    return Scalar(
        data=float(slope),
        title='Flory exponent',
        unit=''
    )
