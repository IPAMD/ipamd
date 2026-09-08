import math
from ipamd.public.constant import na, e, epsilon0, kb
from ipamd.public.utils.output import warning
def func(ff_param, all_info, gala_core):
    rcut = float(ff_param['rcut'])
    epsilon = ff_param['epsilon']
    cs = ff_param['cs']
    temperature = ff_param['temperature']
    neighbor_list = gala_core.NeighborList(all_info, rcut, rcut * 0.1)
    neighbor_list.exclusion(["bond", "angle"])
    if cs < 1e-7:
        cs = 1e-7
    bjerrum_length = e ** 2 / (4 * math.pi * epsilon * epsilon0 * kb * temperature)
    warning(f'Bjerrum length: {bjerrum_length}')
    debye_length = math.sqrt(1 / (8 * math.pi * bjerrum_length * cs * na * 1000)) * 1e9
    debye_force = gala_core.DebyeForce(all_info, neighbor_list, rcut)
    debye_force.setDebyeLength(debye_length)
    if ff_param['shift'] == 'True':
        debye_force.set_energy_shift(True)
    return debye_force
