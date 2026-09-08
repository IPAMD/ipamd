"""update the force field of the molecule"""
configure = {
    "resource": ['ff'],
}
def func(molecule, ff=None):
    """plugin main function"""
    for atom in molecule.atoms:
        atom['prototype'].update(ff)
