"""
plugin for calculating mean squared displacement along a trajectory
"""
import numpy as np
from ipamd.public.models.data import Vector
from ipamd.public.utils.output import warning
from ipamd.public.utils.parser import frame_index_list, protein_range_split

def _trajectory(box, indices, target_molecule, n, per_particle):
    series = []
    for index in indices:
        prop = box.frame(index).current_frame().properties(ignoring_image=True)
        molecules = prop['molecules']
        if n is not None:
            molecules = [molecules[n]]
        else:
            target_molecule_type, target_range = protein_range_split(target_molecule)
            selected = []
            for molecule in molecules:
                if target_molecule_type != "" and molecule['molecule_type'] != target_molecule_type:
                    continue
                if target_range:
                    molecule = {
                        'mass': [molecule['mass'][i] for i in target_range],
                        'position': [molecule['position'][i] for i in target_range]
                    }
                selected.append(molecule)
            molecules = selected

        frame_values = []
        for molecule in molecules:
            positions = np.array(molecule['position'], dtype=float)
            masses = np.array(molecule['mass'], dtype=float)
            if per_particle:
                frame_values.extend(positions)
            else:
                total_mass = np.sum(masses)
                frame_values.append(np.sum(positions * masses[:, None], axis=0) / total_mass)
        series.append(np.array(frame_values, dtype=float))
    return series

def func(box, target_frame, target_molecule=None, n=None, dt=None, per_particle=False, **kwargs):
    indices = frame_index_list(target_frame)
    if len(indices) < 2:
        warning('MSD needs at least two frames')
        return Vector(
            data=[0.0],
            title='MSD',
            x_label='Lag (frames)' if dt is None else 'Time',
            y_label='MSD (nm^2)',
            x_axis=[0]
        )

    series = _trajectory(box, indices, target_molecule, n, per_particle)
    n_frames = len(series)
    msd = np.zeros(n_frames, dtype=float)
    for lag in range(n_frames):
        displacements = []
        for t in range(n_frames - lag):
            delta = series[t + lag] - series[t]
            displacements.append(np.mean(np.sum(delta ** 2, axis=1)))
        msd[lag] = float(np.mean(displacements))

    if dt is None:
        x_axis = list(range(n_frames))
        x_label = 'Lag (frames)'
    else:
        x_axis = [lag * dt for lag in range(n_frames)]
        x_label = 'Time'

    return Vector(
        data=msd,
        title='MSD',
        x_label=x_label,
        y_label='MSD (nm^2)',
        x_axis=x_axis
    )
