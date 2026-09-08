import copy
import os
import sys
from importlib import import_module
import inspect

class PluginBase:
    @staticmethod
    def call(func_name, *args, **kwargs):
        stack = inspect.stack()
        caller_frame = stack[2].frame
        caller_function_object = caller_frame.f_locals.get('self', None)
        func = getattr(caller_function_object, func_name)
        return func(*args, **kwargs)

    def __init__(self, plugin_dir_list):
        self.__schema = {
            'direct': {
            },
            "full": {
                'slf': lambda: self
            }
        }
        self.__resource = {}
        self.__available_plugins = {}
        for plugin_dir in plugin_dir_list:
            available_plugin_file_list = os.listdir(plugin_dir)
            for file in available_plugin_file_list:
                if os.path.isdir(os.path.join(plugin_dir, file)):
                    plugin_path = os.path.join(plugin_dir, file)
                    files = os.listdir(plugin_path)
                    for f in files:
                        if f == 'main.py':
                            self.__available_plugins[file] = {
                                'type': 'dir',
                                'location': plugin_dir,
                                'loaded': False
                            }
                else:
                    filename, extension_name = os.path.splitext(file)
                    if extension_name == '.py':
                        self.__available_plugins[filename] = {
                            'type': 'single_file',
                            'location': plugin_dir,
                            'loaded': False
                        }

    def def_schema(self, name, schema):
        self.__schema[name] = schema

    def add_resource(self, name, value):
        self.__resource[name] = value

    def load_all(self):
        for plugin in self.__available_plugins.keys():
            self.load(plugin)

    def load(self, plugin_name):
        plugin_info = self.__available_plugins[plugin_name]
        sys.path.append(plugin_info['location'])
        if plugin_info['loaded']:
            return
        if plugin_info['type'] == 'dir':
            import_path = plugin_name + '.main'
        else:
            import_path = plugin_name
        module = import_module(import_path)
        plugin_default_config = {
            'schema': 'direct',
            'apply': [],
            'attr': {}
        }
        if hasattr(module, 'configure'):
            configure = module.configure
            for key in plugin_default_config.keys():
                if key not in configure.keys():
                    configure[key] = plugin_default_config[key]
        else:
            configure = plugin_default_config
        schema_name = configure['schema']
        applied_recourses = configure['apply']
        extra_param_dict = copy.deepcopy(self.__schema[schema_name])
        for resource in applied_recourses:
            extra_param_dict[resource] = self.__resource[resource]

        def wrapped_function(self, *args, **kwargs):
            p = {}
            for k in extra_param_dict.keys():
                value = extra_param_dict[k]
                if callable(value):
                    p[k] = value()
                else:
                    p[k] = value
            return module.func(*args, **kwargs, **p)

        wrapped_function.attr = configure['attr']
        setattr(self, plugin_name, wrapped_function.__get__(self, self.__class__))

        plugin_info['loaded'] = True
        sys.path.remove(plugin_info['location'])