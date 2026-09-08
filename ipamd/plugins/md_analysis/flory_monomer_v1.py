"""
plugin for estimating the Flory scaling exponent of a single-chain polymer
"""
import random
import math
import numpy as np
from ipamd.public.models.data import Scalar
from ipamd.public.utils.output import error, warning
from ipamd.public.utils.parser import frame_index_list, protein_range_split

def _calculate_rg2(molecule, start_index, end_index):
    total_mass = 0.0
    mass_center = np.zeros(3, dtype=float)
    for i in range(start_index, end_index + 1):
        total_mass += molecule['mass'][i]
        mass_center += molecule['mass'][i] * np.array(molecule['position'][i], dtype=float)
    mass_center /= total_mass

    r2_m = 0.0
    for i in range(start_index, end_index + 1):
        distance2 = np.linalg.norm(np.array(molecule['position'][i], dtype=float) - mass_center) ** 2
        r2_m += distance2 * molecule['mass'][i]
    return r2_m / total_mass

def func(box, target_frame, target_molecule=None, type_='', sub=1, **kwargs):
    indices = frame_index_list(target_frame)
    if len(indices) == 0:
        indices = [box.current_frame().no]
    if len(indices) <= 5:
        warning("Small number of frames may lead to inaccurate results.")

    selector = target_molecule if target_molecule not in (None, '') else type_
    target_molecule_type, target_range = protein_range_split(selector)

    log_n_list = []
    log_r_list = []
    rg_full_molecule = []
    molecule_length = 0

    for index in indices:
        prop = box.frame(index).current_frame().properties(ignoring_image=True)
        molecules = prop['molecules']
        if target_molecule_type == '' and molecules:
            target_molecule_type = molecules[0]['molecule_type']
        for molecule in molecules:
            if molecule['molecule_type'] != target_molecule_type:
                continue
            if target_range:
                molecule = {
                    'mass': [molecule['mass'][i] for i in target_range],
                    'position': [molecule['position'][i] for i in target_range],
                    'n_atoms': len(target_range)
                }
            molecule_length = molecule['n_atoms']
            if molecule_length < 30:
                error('Molecule too small')
                return Scalar(data=-1, title='Flory exponent', unit='')
            sampling_times = int(math.pow(math.log(molecule_length), 2) * sub)
            for _ in range(sampling_times):
                nums = np.arange(10, molecule_length + 1)
                weight = np.linspace(10, molecule_length + 1, molecule_length - 9)
                weight = np.divide(1, weight)
                target_length = np.random.choice(nums, p=weight / np.sum(weight))
                start_index = random.randint(0, molecule_length - target_length)
                end_index = start_index + target_length - 1
                rg2 = _calculate_rg2(molecule, start_index, end_index)
                log_n_list.append(math.log(target_length))
                log_r_list.append(math.log(rg2) / 2)
            rg_full_molecule.append(_calculate_rg2(molecule, 0, molecule_length - 1))

    if len(log_n_list) == 0:
        error('No molecule matched the selected type')
        return Scalar(data=0, title='Flory exponent', unit='')

    sorted_index = np.argsort(log_n_list)
    log_n_list = np.array(log_n_list)[sorted_index]
    log_r_list = np.array(log_r_list)[sorted_index]
    max_log_length = math.log(molecule_length)
    split_index = np.searchsorted(log_n_list, max_log_length * 2 / 3, side='right')
    sub_n_list = log_n_list[:int(split_index)]
    sub_r_list = log_r_list[:int(split_index)]
    design_short = np.vstack([sub_n_list, np.ones(len(sub_n_list))]).T
    _, intercept_short = np.linalg.lstsq(design_short, sub_r_list, rcond=None)[0]
    design_all = np.vstack([log_n_list, np.ones(len(log_n_list))]).T
    slope_all, _ = np.linalg.lstsq(design_all, log_r_list, rcond=None)[0]
    if slope_all > 0.53:
        exponent = slope_all
    else:
        rg_molecule = np.average(rg_full_molecule)
        exponent = (math.log(rg_molecule) / 2 - intercept_short) / math.log(molecule_length)

    return Scalar(
        data=float(exponent),
        title='Flory exponent',
        unit=''
    )
