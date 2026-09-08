import os
import pypdbio
from ipamd.public.utils.output import error, info
configure = {
    "resource": ['persistency_dir'],
}
def func(pdb_id, persistency_dir=None):
    try:
        pypdbio.fetch(
            pdb_id,
            os.path.join(persistency_dir, pdb_id + '.pdb'),
        )
    except Exception as e:
        error('Failed to download ' + pdb_id + ' from RCSB. ' + str(e))
        return
    info('Downloaded ' + pdb_id + ' from RCSB.')
