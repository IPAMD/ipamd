from ipamd.public.models.md import Box
from ipamd.public.utils.output import warning
configure = {
    "resource": ['persistency_dir', 'ff'],
}

def func(filename, ff, persistency_dir=None):
    tmp_box = Box(0, 0, 0, ff, persistency_dir)
    tmp_box.new_frame()
    tmp_box.read_xml(filename)
    molecules = tmp_box.current_frame().molecules
    if len(molecules) > 1:
        warning('More than 1 molecule detected, only the first one will be parsed.')
    tmp_box.current_frame().molecules = [molecules[0]]
    tmp_box.align_center()
    return tmp_box.current_frame().molecules[0]['prototype']
