"""
plugin to calculate the free energy of the slab simulation
"""
import math
import copy
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
import numpy as np
from ipamd.public.utils.plugin_manager_v1 import PluginBase
from ipamd.public.constant import r
from ipamd.public.models.data import Vector
from ipamd.public.utils.output import warning, error


def func(box, target_frame, direction='Z', d=1, window=5, **kwargs):
    """
    Calculate the free energy of the slab simulation
    :param box: the box of the simulation, don't need to be set
    :param target_frame: the target frame of the simulation, the frame to calculate the free energy
    :param direction: the direction of the simulation
    :param d: divide the simulation into parts with width d
    :param window: sliding window size for density smoothing; 1 disables smoothing
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

    data = copy.deepcopy(res.data)
    if window > 1:
        data = uniform_filter1d(np.asarray(data, dtype=float), size=window, mode='wrap')

    max_density = data.max()
    min_density = data.min()

    first_cross_index = None
    last_cross_index = None
    for i in range(len(data)-1):
        if data[i] < min_density + (max_density - min_density) * 0.75 <= data[i + 1]:
            if first_cross_index is None:
                first_cross_index = i
        elif data[i] >= min_density + (max_density - min_density) * 0.75 > data[i + 1]:
            last_cross_index = i + 1


    center = (first_cross_index + last_cross_index) / 2
    data = res.data

    left_side_data = data[:int(center) + 1][::-1]
    right_side_data = data[math.ceil(center):]

    def phase_boundary_func(x, rho1, rho2, d, z0):
        return 0.5 * (rho1 + rho2) - 0.5 * (rho1 - rho2) * np.tanh((x - z0) / d)

    def generate_x_list(data_length, d):
        return [i * d + 0.5 * d for i in range(data_length)]

    def fit_phase_boundary(y):
        rho1_guess = y.max()
        rho2_guess = y.min()
        d_guess = d
        z0_guess = d * len(y) / 2
        rho1_fit, rho2_fit, d_fit, z0_fit = curve_fit(
            phase_boundary_func,
            np.array(generate_x_list(len(y), d)),
            y,
            p0=[rho1_guess, rho2_guess, d_guess, z0_guess]
        )[0]
        if d_fit >= 0:
            return rho1_fit, rho2_fit, d_fit, z0_fit
        else:
            return rho2_fit, rho1_fit, -d_fit, z0_fit
    rho11, rho12, d1, z01 = fit_phase_boundary(left_side_data)
    rho21, rho22, d2, z02 = fit_phase_boundary(right_side_data)
    previous_z01 = 0
    previous_z02 = 0

    n_iter = 0
    max_iter = 20
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
        rho11, rho12, d1, z01 = fit_phase_boundary(extended_left_side_data)
        extended_right_side_data = np.concatenate((
            extended_right_side_data_condense_phase,
            right_side_data,
            extended_right_side_data_dilute_phase
        ))
        rho21, rho22, d2, z02 = fit_phase_boundary(extended_right_side_data)
        z01 = z01 - d * len(extended_left_side_data_condense_phase)
        z02 = z02 - d * len(extended_right_side_data_condense_phase)
        n_iter += 1
        if n_iter > max_iter:
            error('The free energy calculation failed after 20 iterations. Please check the data.')
            return None

    rho_bar1 = 0.5 * (rho11 + rho21)
    if rho22 < 1e-8:
        rho22 = 1e-8
    if rho12 < 1e-8:
        rho12 = 1e-8
    rho_bar2 = 0.5 * (rho12 + rho22)


    if 0.67 < rho_bar1 / rho_bar2 < 1.5:
        warning(
            'The difference of density between the two phases is too small. '
            'Free energy calculation may be inaccurate.'
        )

    temperature = box.env.values['temperature']
    free_energy = - r * temperature * math.log(rho_bar1 / rho_bar2) / 1000
    return Vector(
        title='Free Energy',
        x_axis=['Free Energy (kJ/mol)', 'density(dense) (g/mL)', 'density(dilute) (g/mL)'],
        data=[free_energy, rho_bar1, rho_bar2]
    )
