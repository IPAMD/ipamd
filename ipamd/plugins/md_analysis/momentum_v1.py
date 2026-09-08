"""
plugin for calculating the total momentum of the system
"""
from ipamd.public.constant import na
from ipamd.public.models.data import Vector

configure = {
    "schema": ['frame'],
}

def func(frame, **kwargs):
    prop = frame.properties(ignoring_image=False)
    molecules = prop['molecules']
    px = 0.0
    py = 0.0
    pz = 0.0
    for molecule in molecules:
        mass_list = molecule['mass']
        vel_list = molecule['velocity']
        n_atoms = len(mass_list)
        for i in range(n_atoms):
            mass_kg = mass_list[i] / na / 1000
            velocity = vel_list[i] if vel_list[i] is not None else (0.0, 0.0, 0.0)
            velocity_m_s = [v * 1000 for v in velocity]
            px += mass_kg * velocity_m_s[0]
            py += mass_kg * velocity_m_s[1]
            pz += mass_kg * velocity_m_s[2]

    return Vector(
        data=[px, py, pz],
        title='Momentum',
        x_label='Component',
        y_label='Momentum (kg m/s)',
        x_axis=['px', 'py', 'pz']
    )
