def func(ff_param, all_info, gala_core):
    rcut = float(ff_param['rcut'])
    neighbor_list = gala_core.NeighborList(all_info, rcut, rcut * 0.2)
    neighbor_list.exclusion(["bond", "angle"])
    ah_force = gala_core.AHForce(all_info, neighbor_list, rcut)
    if ff_param['shift'] == 'True':
        ah_force.set_energy_shift(True)
    atom_type = all_info.getBasicInfo().getParticleTypes()
    epsilon = 0.8368
    num_atom_types = len(atom_type)
    for i in range(num_atom_types):
        for j in range(i, num_atom_types):
            atom1 = atom_type[i]
            atom2 = atom_type[j]
            # atom not defined
            if atom1 not in ff_param or atom2 not in ff_param:
                continue

            lambda_ = 0.5 * (float(ff_param[atom1]['lambda']) + float(ff_param[atom2]['lambda']))
            sigma = 0.5 * (float(ff_param[atom1]['sigma']) + float(ff_param[atom2]['sigma']))
            ah_force.setParams(atom1, atom2, epsilon, sigma, lambda_)
    return ah_force
