def func(ff_param, all_info, gala_core):
    stacking_force = gala_core.AHPairForce(all_info)
    all_bond_types = all_info.getBondInfo().getBondTypes()
    stacking_force.setRcut(float(ff_param['rcut']))
    stacking_force.setMultiplyFactor(15)
    stacking_force.setEnergyShift(ff_param['shift'] == 'True')
    if 'RS-RS' in all_bond_types:
        stacking_force.setParams(
            "RS-RS",
            0.8368,
            1.18,
            0.4
        )
    return stacking_force
