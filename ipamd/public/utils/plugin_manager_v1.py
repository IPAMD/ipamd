import os
from importlib import util as importlib_util
from ipamd.public.utils.output import tabulate
import inspect

class PluginBase:
    @staticmethod
    def call(func_name, *args, module=None, **kwargs):
        stacks = inspect.stack()
        caller_function_object = None
        for stack in stacks:
            if stack.function == 'ipamd_wrapped_function':
                caller_function_object = stack.frame.f_locals.get('self', None)
                if module is not None:
                    caller_function_object = getattr(caller_function_object.app, module)
                break
        func = getattr(caller_function_object, func_name)
        return func(*args, **kwargs)

    def __init__(self, plugin_dir_list):
        self.__schema = {}
        self.__resource = {}
        self.__available_plugins = {}
        for plugin_dir in plugin_dir_list:
            available_plugin_file_list = os.listdir(plugin_dir)
            for file in available_plugin_file_list:
                filename, extension_name = os.path.splitext(file)
                if extension_name == '.py':
                    self.__available_plugins[filename] = {
                        'location': plugin_dir,
                        'loaded': False
                    }

    def def_schema(self, name, schema):
        self.__schema[name] = schema

    def add_resource(self, name, value):
        self.__resource[name] = value

    def plugin_info(self):
        tabulate(
            title="Plugin Information",
            headers=['Name','Location', 'Status'],
            rows=[
                (plugin_name,
                 self.__available_plugins[plugin_name]['location'],
                 'loaded' if self.__available_plugins[plugin_name]['loaded'] else '')
                for plugin_name in self.__available_plugins
            ]
        )

    def load_all(self):
        for plugin in self.__available_plugins:
            self.load(plugin)

    def load(self, plugin_name):
        plugin_info = self.__available_plugins[plugin_name]
        if plugin_info['loaded']:
            return

        plugin_path = os.path.join(plugin_info['location'], plugin_name + '.py')
        unique_name = f"ipamd_{os.path.basename(plugin_info['location'])}_{plugin_name}"
        spec = importlib_util.spec_from_file_location(unique_name, plugin_path)
        module = importlib_util.module_from_spec(spec)
        spec.loader.exec_module(module)
        configure = {
            "type": "function",
            "schema": [],
            "resource": [],
            "alias": None,
            "attributes": {}
        }

        if hasattr(module, 'configure'):
            for key in module.configure:
                configure[key] = module.configure[key]

        schema = configure['schema']
        recourses = configure['resource']

        if configure['type'] == 'function':
            def ipamd_wrapped_function(self, *args, **kwargs):
                p = {}
                for s in schema:
                    p.update(self.__schema[s](*args, **kwargs))
                for res in recourses:
                    if res in self.__resource:
                        if callable(self.__resource[res]):
                            p[res] = self.__resource[res]()
                        else:
                            p[res] = self.__resource[res]
                return module.func(*args, **kwargs, **p)

            ipamd_wrapped_function.attr = configure['attributes']
            function_name = configure['alias'] if configure['alias'] else plugin_name
            setattr(self, function_name, ipamd_wrapped_function.__get__(self, self.__class__))

        elif configure['type'] == 'class':
            class WrappedClass(module.Cls):
                def __init__(self, *args, **kwargs):
                    super().__init__(*args, **kwargs)
                    for res in recourses:
                        if res in self.__resource:
                            setattr(self, res, self.__resource[res])
                    for attr in configure['attributes']:
                        setattr(self, attr, configure['attributes'][attr])

            class_name = configure.get('alias', plugin_name)
            setattr(self, class_name, WrappedClass)
        else:
            raise ValueError(f"Unknown plugin type: {configure['type']}")
        plugin_info['loaded'] = True
