from ipamd.public.utils.plugin_manager_v1 import PluginBase
from ipamd.public import shared_data

class MdAnalysis(PluginBase):
    def __init__(self, app):
        plugin_dir = [shared_data.module_installation_dir + '/plugins/md_analysis']
        super().__init__(plugin_dir)
        self.def_schema(
            'frame', 
            lambda **kwargs: {
                'frame': kwargs['box'].frame(kwargs['target_frame']).current_frame()
            }
        )
        self.app = app
        self.add_resource('working_dir', lambda: self.app.working_dir)
        self.load_all()
