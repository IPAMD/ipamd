def func(ff_param, all_info, gala_core):
    bond_force = gala_core.BondForceHarmonic(all_info)
    all_bond_types = all_info.getBondInfo().getBondTypes()
    for bond_type in ff_param.keys():
        if bond_type in all_bond_types:
            bond_force.setParams(bond_type, float(ff_param[bond_type]['k']), float(ff_param[bond_type]['r0']))
    return bond_force
