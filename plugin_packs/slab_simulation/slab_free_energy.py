"""
plugin to calculate the free energy of the slab simulation
"""
import math
import sys
from scipy.optimize import curve_fit
import numpy as np
from ipamd.public.utils.plugin_manager_v1 import PluginBase
from ipamd.public.constant import r
from ipamd.public.models.data import Vector
from ipamd.public.utils.output import warning, error


def func(box, target_frame, direction='Z', d=1, **kwargs):
    """
    Calculate the free energy of the slab simulation
    :param box: the box of the simulation, don't need to be set
    :param target_frame: the target frame of the simulation, the frame to calculate the free energy
    :param direction: the direction of the simulation
    :param d: divide the simulation into parts with width d
    :param kwargs: other parameters, don't need to be set

    :return: a vector containing the free energy, the density of two phases
    """
    res = PluginBase.call(
        'batch_compute',
        'density_align',
        box=box,
        target_frame=target_frame,
        direction=direction,
        d=d
    )
    data = res.data
    max_density = data.max()
    first_cross_index = None
    last_cross_index = None
    for i in range(len(data)-1):
        if data[i] < max_density * 0.5 <= data[i + 1]:
            if first_cross_index is None:
                first_cross_index = i
        elif data[i] >= max_density * 0.5 > data[i + 1]:
            last_cross_index = i + 1

    try:
        center = (first_cross_index + last_cross_index) / 2
    except TypeError:
        error('Falied to find the phase boundary. Please check the simulation settings.')
        sys.exit(1)

    left_side_data = data[:int(center) + 1][::-1]
    right_side_data = data[math.ceil(center):]

    def phase_boundary_func(x, rho1, rho2, d, z0):
        return 0.5 * (rho1 + rho2) - 0.5 * (rho1 - rho2) * np.tanh((x - z0) / d)

    def generate_x_list(data_length, d):
        x_list = []
        for i in range(data_length):
            x_list.append(i * d + 0.5 * d)
        return x_list

    left_x_list = generate_x_list(len(left_side_data), d)

    rho11, rho12, d1, z01 = curve_fit(
        phase_boundary_func,
        np.array(left_x_list),
        left_side_data
    )[0]

    right_x_list = generate_x_list(len(right_side_data), d)
    rho21, rho22, d2, z02 = curve_fit(
        phase_boundary_func,
        np.array(right_x_list),
        right_side_data
    )[0]

    previous_z01 = 0
    previous_z02 = 0

    while abs(z01 - previous_z01) > d or abs(z02 - previous_z02) > d:
        previous_z01 = z01
        previous_z02 = z02
        extended_left_side_data_condense_phase = right_side_data[
            0:max(int((z02 - 2.5 * d2) / d), 0)
        ][::-1]
        extended_left_side_data_dilute_phase = right_side_data[
            min(int((z02 + 2.5 * d2) / d), len(right_side_data)):
        ][::-1]
        extended_left_side_data = np.concatenate((
            extended_left_side_data_condense_phase,
            left_side_data,
            extended_left_side_data_dilute_phase
        ))

        extended_right_side_data_condense_phase = left_side_data[
            0:max(int((z01 - 2.5 * d1) / d), 0)
        ][::-1]
        extended_right_side_data_dilute_phase = left_side_data[
            min(int((z01 + 2.5 * d1) / d), len(left_side_data)):
        ][::-1]
        extended_left_x_list = generate_x_list(len(extended_left_side_data), d)
        extended_right_side_data = np.concatenate((
            extended_right_side_data_condense_phase,
            right_side_data,
            extended_right_side_data_dilute_phase
        ))
        extended_right_x_list = generate_x_list(len(extended_right_side_data), d)
        rho11, rho12, d1, z01 = curve_fit(
            phase_boundary_func,
            np.array(extended_left_x_list),
            extended_left_side_data
        )[0]
        z01 = z01 - d * len(extended_left_side_data_condense_phase)
        rho21, rho22, d2, z02 = curve_fit(
            phase_boundary_func,
            np.array(extended_right_x_list),
            extended_right_side_data
        )[0]
        z02 = z02 - d * len(extended_right_side_data_condense_phase)

    rho_bar1 = 0.5 * (rho11 + rho21)
    if rho22 < 1e-8:
        rho22 = 1e-8
    if rho12 < 1e-8:
        rho12 = 1e-8
    rho_bar2 = 0.5 * (rho12 + rho22)
    if 0.67 < rho_bar1 / rho_bar2 < 1.5:
        warning('Free energy calculation may be inaccurate.')

    temperature = box.env.values['temperature']
    free_energy = - r * temperature * math.log(rho_bar1 / rho_bar2) / 1000
    return Vector(
        title='Free Energy',
        x_axis=['Free Energy (kJ/mol)', 'density(dense) (g/mL)', 'density(dilute) (g/mL)'],
        data=[free_energy, (rho11 + rho21) / 2, (rho12 + rho22) / 2]
    )
