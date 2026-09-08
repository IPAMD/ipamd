"""
plugin for calculating the instantaneous temperature of the system
"""
from ipamd.public.constant import na, kb
from ipamd.public.models.data import Scalar

configure = {
    "schema": ['frame'],
}

def func(frame, **kwargs):
    prop = frame.properties(ignoring_image=False)
    molecules = prop['molecules']
    ek = 0.0
    n_free = 0
    rigid_groups = set()
    for molecule_index, molecule in enumerate(molecules):
        mass_list = molecule['mass']
        vel_list = molecule['velocity']
        rigid_list = molecule['rigid_group']
        n_atoms = len(mass_list)
        for i in range(n_atoms):
            mass_kg = mass_list[i] / na / 1000
            velocity = vel_list[i] if vel_list[i] is not None else (0.0, 0.0, 0.0)
            velocity_m_s = [v * 1000 for v in velocity]
            ek += 0.5 * mass_kg * (
                velocity_m_s[0] ** 2 + velocity_m_s[1] ** 2 + velocity_m_s[2] ** 2
            )
            rigid = rigid_list[i]
            if rigid == -1:
                n_free += 1
            else:
                rigid_groups.add((molecule_index, rigid))

    freedom = 3 * n_free + 6 * len(rigid_groups)
    temperature = 2 * ek / (freedom * kb) if freedom > 0 else 0.0
    return Scalar(
        data=temperature,
        title='Temperature',
        unit='K'
    )
