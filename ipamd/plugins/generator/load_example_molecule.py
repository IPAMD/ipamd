import os
from shutil import copyfile
from ipamd.public.utils.output import error
from ipamd.public import shared_data
configure = {
    "resource": ['persistency_dir']
}

def func(molecule_name, persistency_dir=None):
    mol_db_dir = os.path.join(shared_data.module_installation_dir, 'data/molecules')
    molecule_file_path = os.path.join(mol_db_dir, molecule_name + '.pdb')
    target_file_path = os.path.join(persistency_dir, molecule_name + '.pdb')
    if os.path.exists(molecule_file_path):
        copyfile(
            molecule_file_path,
            target_file_path
        )
    else:
        error(f"No molecule {molecule_name} file found.")
        return
