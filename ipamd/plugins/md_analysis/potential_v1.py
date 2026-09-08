"""
plugin for reading potential energy from a simulation log
"""
import os
from ipamd.public.models.data import Scalar
from ipamd.public.utils.output import error

configure = {
    "schema": ['frame'],
    "resource": ['working_dir'],
}

def func(frame, working_dir, simulation, **kwargs):
    frame_no = frame.no
    simulation_name = simulation.job_name
    log_file = os.path.join(working_dir, simulation_name + '.log')
    if not os.path.exists(log_file):
        error('simulation should be run first')
        raise FileNotFoundError(log_file)

    time_step = frame_no * simulation.period
    potential = None
    with open(log_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    for line in lines[1:]:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split()
        time_step_of_line = float(parts[0])
        if abs(time_step_of_line - time_step) < 1e-6:
            potential = float(parts[3])
            break

    if potential is None:
        error('simulation should be run first')
        raise ValueError(f'potential energy at timestep {time_step} not found in {log_file}')

    return Scalar(
        data=potential,
        title='Potential Energy',
        unit='kJ/mol'
    )
