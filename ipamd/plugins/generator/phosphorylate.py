import copy
from ipamd.public.utils.output import warning

configure = {
    "resource": ['ff'],
}
def func(molecule, index, ff=None):
    mol = copy.deepcopy(molecule)
    ff = copy.deepcopy(ff)

    target_atom = mol.atoms[index]
    pt = target_atom['prototype']
    typ = pt.get('type')
    if typ != 'S' and typ != 'T' and typ != 'Y':
        warning('Only SER, THR, and TYR can be phosphorylated.')

    atom_definition = ff.atom_definition[typ]
    if atom_definition['charge'].startswith('compute:'):
        atom_definition['charge'] += '-1.0'
    else:
        atom_definition['charge'] = str(float(atom_definition['charge']) - 1.0)

    if atom_definition['mass'].startswith('compute:'):
        atom_definition['mass'] += '+79.982'
    else:
        atom_definition['mass'] = str(float(atom_definition['mass']) + 79.982)

    ah_parameter = ff.ff_param['ah'][typ]

    extra_ff = {
        "atom_definition": {
            f"{typ}P": atom_definition
        },
        "ff_param": {
            "ah": {
                f"{typ}P": ah_parameter
            }
        }
    }

    actual_charge = target_atom.get('charge')
    actual_mass = target_atom.get('mass')
    if actual_mass is not None:
        if actual_mass.startswith('compute:'):
            actual_mass += '+79.982'
        else:
            actual_mass = str(float(actual_mass) + 79.982)
        target_atom.set('mass', actual_mass)
    if actual_charge is not None:
        if actual_charge.startswith('compute:'):
            actual_charge += '-1.0'
        else:
            actual_charge = str(float(actual_charge) - 1.0)
        target_atom.set('charge', actual_charge)

    pt.set('type', f'{typ}P')

    return mol, extra_ff