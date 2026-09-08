import random
def func(param, all_info, group, gala_core):
    t_reduced = param['temperature'] * 8.3143 / 1000.0
    rigid = param['rigid']
    integrator = None
    if rigid:
        integrator = gala_core.LangevinNVTRigid(all_info, group, t_reduced, random.randint(1, 100))
    else:
        integrator = gala_core.LangevinNVT(all_info, group, t_reduced, random.randint(1, 100))
    integrator.setGamma(0.01)
    return integrator
