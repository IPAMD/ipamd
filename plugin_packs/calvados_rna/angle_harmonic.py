def func(ff_param, all_info, gala_core):
    angle_force = gala_core.AngleForceHarmonic(all_info)
    all_angle_types = all_info.getAngleInfo().getAngleTypes()
    for angle_type in ff_param.keys():
        if angle_type in all_angle_types:
            angle_force.setParams(
                angle_type, float(ff_param[angle_type]['k']), float(ff_param[angle_type]['theta0'])
            )
    return angle_force
