from ipamd.public.models.md import Box
configure = {
    "resource": ['persistency_dir', 'ff'],
}

def func(filename, ff, persistency_dir=None):
    box = Box(0, 0, 0, ff, persistency_dir)
    box.new_frame()
    box.read_xml(filename) # pylint: disable=no-member
    return box
