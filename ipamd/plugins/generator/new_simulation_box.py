from ipamd.public.models.md import Box
configure = {
    "resource": ['ff', 'persistency_dir'],
}
def func(x, y, z, ff=None, persistency_dir=None):
    return Box(x, y, z, ff, persistency_dir).new_frame()
