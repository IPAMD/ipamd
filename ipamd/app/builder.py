from ipamd.public.utils.plugin_manager_v1 import PluginBase
from ipamd.public import shared_data

class Builder(PluginBase):
    def __init__(self, app):
        extra_plugin_dir = shared_data.config.get('builder_plugin_dir')
        plugin_dir = [shared_data.module_installation_dir + '/plugins/generator'] + extra_plugin_dir
        super().__init__(plugin_dir)
        self.__app = app
        self.add_resource('persistency_dir', lambda: self.__app.working_dir)
        self.add_resource('cuda', self.__app.cuda)
        self.add_resource('ff', lambda: self.__app.force_field)
        if shared_data.config.get('auto_load'):
            self.load_all()
