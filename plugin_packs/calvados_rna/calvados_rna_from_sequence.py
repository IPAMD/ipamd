import numpy as np
from ipamd.public.utils.plugin_manager_v1 import PluginBase
from ipamd.public.models.sequence import RNASequence
from ipamd.public.utils.output import warning
from ipamd.public.models.md import Atom

configure = {
    "resource": ['ff'],
}

def func(mol: RNASequence, curve_func, reverse=False, fallback_ref_vec=[1, 1, 1], ff=None):
    """
    :param mol: the RNA sequence object
    :param curve_func: the curve function
    :param fallback_ref_vec: the residue will be placed on the same plane as the fallback_ref_vec, only when all bonds are parallel this parameter is used.
    :return: the molecule object
    """
    molecule = PluginBase.call(
        "molecule_from_curve",
        mol,
        curve_func,
        start_add=17.007,
        end_add=1.008,
        new_bond="RB-RB",
        rename_map={
            'U': 'NMPB',
            'A': 'NMPB',
            'C': 'NMPB',
            'G': 'NMPB',
        },
    )
    bond_vec_list = []

    last_atom = None
    for atom in molecule.atoms:
        if last_atom is not None:
            bond_vec_list.append(np.array(atom['offset']) - np.array(last_atom['offset']))
        last_atom = atom

    # check if all bonds are parallel
    all_parallel = True
    ref_bond_vec = bond_vec_list[0]
    for bond_vec in bond_vec_list:
        cross_product = np.cross(ref_bond_vec, bond_vec)
        if np.linalg.norm(cross_product) > 1e-6:
            all_parallel = False
            break
    if all_parallel:
        warning('All bonds are parallel.')

    # calculate the position of side chain
    bond_vec_list = [-bond_vec_list[0]] + bond_vec_list + [-bond_vec_list[-1]]

    def normalize(vec):
        return vec / np.linalg.norm(vec)
    for index, atom in enumerate(molecule.atoms):
        if index == len(mol):
            break

        vec1 = normalize(bond_vec_list[index])
        vec2 = normalize(bond_vec_list[index + 1])
        ref_vec = vec2
        if all_parallel:
            ref_vec = fallback_ref_vec
        else:
            if np.linalg.norm(np.cross(vec1, vec2)) < 1e-6:
                for vec in bond_vec_list:
                    if np.linalg.norm(np.cross(vec, vec2)) > 1e-6:
                        ref_vec = vec
                        break

        # make sure the angle between vec1 and vec2 is greater than 90 degrees
        dot_product = np.dot(vec1, vec2)
        if dot_product > 0:
            vec1 = -vec1

        if np.linalg.norm(np.cross(vec1, vec2)) > 1e-6:
            # vec1 and vec2 are not parallel
            target_vec = normalize(vec1 + vec2)
        else:
            # vec1 and vec2 are parallel
            target_vec = normalize(np.cross(vec1, np.cross(vec1, ref_vec)))
        bond_length = 0.45
        sign = -1 if reverse else 1
        new_offset = [
            atom['offset'][0] + target_vec[0] * bond_length * sign,
            atom['offset'][1] + target_vec[1] * bond_length * sign,
            atom['offset'][2] + target_vec[2] * bond_length * sign,
        ]
        molecule.add_atom(
            atom=Atom(
                velocity=[0, 0, 0],
                atom_type=f'NMP{mol[index].upper()}',
                ff=ff
            ),
            coordinate=new_offset,
            bond=None
        )
        molecule.link(index, index + len(mol), 'RB-RS')
        if index > 0:
            molecule.link(
                index + len(mol), index + len(mol) - 1, 'RS-RS'
            )
        if index > 0 and index < len(mol) - 1:
            molecule.add_angle(index - 1, index, index + 1, 'RB-RB-RB')
    return molecule
