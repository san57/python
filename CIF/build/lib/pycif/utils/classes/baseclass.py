from __future__ import print_function

import importlib
from inspect import ismethod
from types import MethodType, FunctionType, ModuleType

import numpy as np
import pandas as pd
import yaml
from builtins import object
from builtins import zip

from pycif.utils.check.errclass import PluginError


class Plugin(object):
    """Plugin class

        It is used to store all information about sub-parts of PYVAR
        and loads if needed sub-modules
    """

    # Authorized Plugins
    plugin_types = {
        "mode": [".modes", "Mode"],
        "meteo": [".meteos", "Meteo"],
        "model": [".models", "Model"],
        "fields": [".fields", "Fields"],
        "fluxes": [".fluxes", "Fluxes"],
        "domain": [".domains", "Domain"],
        "obsvect": [".obsvects", "ObsVect"],
        "platform": [".platforms", "Platform"],
        "statevect": [".statevects", "StateVect"],
        "minimizer": [".minimizers", "Minimizer"],
        "simulator": [".simulators", "Simulator"],
        "obsparser": [".obsparsers", "ObsParser"],
        "transform": [".transforms", "Transform"],
        "chemistry": [".chemistries", "Chemistry"],
        "obsoperator": [".obsoperators", "ObsOperator"],
        "measurements": [".measurements", "Measurement"],
    }

    # Registered and already loaded plugins
    registered = {}
    loaded_instances = {}
    reference_instances = {}
    subreference_instances = {}

    def __init__(self, plg_orig=None):
        """Initializes the Plugin.

        Args:
            plg_orig (Plugin): if specified, initialize a new plugin from
            attributes of the origin plugin

        """

        attributes = self._get_attribute_list(plg_orig)

        for attr in attributes:
            setattr(self, attr, getattr(plg_orig, attr))

        if not hasattr(self, "attributes"):
            self.attributes = []

        # Updating requirements in any
        self.default_requirements = {
            key: {"any": True, "empty": True}
            for key in ["datei", "datef", "workdir", "logfile", "verbose"]
        }

        if not hasattr(self, "requirements"):
            self.requirements = {}

    @classmethod
    def _get_attribute_list(cls, plg):
        """Get the list of attributes excluding special methods"""
        return [a for a in dir(plg) if not a.startswith("_")]

    @staticmethod
    def plugin_key(name, version, plugin_type):
        """Creates the key for a name, version and plugin type

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin

        Returns:
            frozenset: dict key for a plugin and version

        Notes:
            This function uses a frozenset to protect the plugin name from any
            later alteration.
        """
        return name, version, plugin_type

    @classmethod
    def print_registered(
        cls, print_requirement=False, types=[], names=[], versions=[]
    ):
        """Print in a user-friendly format the list of available plugins"""

        keys = [
            k
            for k in list(cls.registered.keys())
            if (k[0] in names or names == [])
            and (k[1] in versions or versions == [])
            and (k[2] in types or types == [])
        ]
        names, versions, types = list(zip(*keys))

        modules = [
            cls.get_registered(n, v, t)
            for n, v, t in zip(names, versions, types)
        ]

        print(
            "List of all available plugins and requirements "
            "for each class: \n\n"
        )
        for tt in set(types):
            print("\t", tt)
            for n, v, t, mod in zip(names, versions, types, modules):
                if t == tt:
                    print("\t\t- {}, {}".format(n, v))

                    # Print requirements
                    if print_requirement and hasattr(mod, "requirements"):
                        print("\t\t\tRequires:")
                        for req in mod.requirements:
                            name = mod.requirements[req].get("name", "")
                            version = mod.requirements[req].get("version", "")
                            any = mod.requirements[req].get("any", False)
                            empty = mod.requirements[req].get("empty", "")
                            print("\t\t\t\t- {}:".format(req))
                            print("\t\t\t\t\t- name: {}".format(name))
                            print("\t\t\t\t\t- version: {}".format(version))
                            print("\t\t\t\t\t- any: {}".format(any))
                            print("\t\t\t\t\t- empty: {}".format(empty))

            print("\n")

    @classmethod
    def is_registered(cls, name, version, plugin_type):
        """Check if a plugin is registered

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin

        Returns:
            bool: True if a parser is registered
        """
        return cls.plugin_key(name, version, plugin_type) in cls.registered

    @classmethod
    def is_allowed(cls, plugin_type):
        """Check whether a plugin type is allowed in pyCIF or not

        Args:
            plugin_type (str): type of plugin

        Returns:
            bool: True if allowed plugin
        """

        return plugin_type in cls.plugin_types

    @classmethod
    def is_loaded(cls, name, version, plugin_type):
        """Check whether a plugin is loaded

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin

        Returns:
            bool
        """

        return (
            cls.plugin_key(name, version, plugin_type) in cls.loaded_instances
        )

    @classmethod
    def register_plugin(cls, name, version, module, plugin_type="", **kwargs):
        """Register a module for a plugin and version with possibly options

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin
            module (types.ModuleType): module defining the interface
                between pyCIF and the plugin
            **kwargs (dictionary): default options for module

        """
        if cls.is_registered(name, version, plugin_type):
            raise ValueError(
                "Already created a Module "
                "for plugin {} (version {}) and type {}".format(
                    name, version, plugin_type
                )
            )

        cls.registered[
            cls.plugin_key(name, version, plugin_type)
        ] = module.__name__

    @classmethod
    def get_registered(cls, name, version, plugin_type):
        """Get the correct registered plugin, given its name, version and type

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin

        Returns:
            Plugin: plugin module for plugin version
        """
        if not cls.is_registered(name, version, plugin_type):
            raise PluginError(
                "No {} module found for plugin {} and Version {}".format(
                    plugin_type, name, version
                )
            )

        module = importlib.import_module(
            cls.registered[cls.plugin_key(name, version, plugin_type)]
        )

        return module

    @classmethod
    def get_loaded(cls, name, version, plugin_type):
        """Get the correct loaded plugin, given its name, version and type

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            plugin_type (str): type of plugin

        Returns:
            Plugin: plugin module for plugin version
        """
        if not cls.is_loaded(name, version, plugin_type):
            raise PluginError(
                "No {} module loaded for plugin {} and Version {}".format(
                    plugin_type, name, version
                )
            )

        return cls.loaded_instances[cls.plugin_key(name, version, plugin_type)]

    @classmethod
    def get_subclass(cls, plg_type):
        """Get the plugin class template from a given type

        Args:
            plg_type (str): the plugin type to load
        """

        subclass = cls.plugin_types[plg_type]
        return getattr(
            importlib.import_module(subclass[0], "pycif.utils.classes"),
            subclass[1],
        )

    @classmethod
    def load_registered(cls, name, version, plg_type, plg_orig=None):
        """Get a sub-class instance of a registered plugin.
        This can be used to get default required plugins


        """

        plgtmp = cls.get_registered(name, version, plg_type)

        if plg_orig is not None:
            for attr in plg_orig.attributes:
                setattr(plgtmp, attr, getattr(plg_orig, attr))

            plgtmp.attributes = plg_orig.attributes[:]

        # Adding plugin attribute
        plgtmp.plugin = cls.from_dict(
            {"name": name, "version": version, "type": plg_type}
        )

        # Creating a sub-class instance and initializing it
        plgtmp = cls.childclass_factory(plg_orig=plgtmp)

        # Initializing the sub-class instance
        plgtmp.initiate_template()

        return plgtmp

    def _load_plugin_type(self, key, parent_plg_type=None):
        """Load plugin type and add it to the plugin if not already specified

        Args:
            key (str): the plugin type
            parent_plg_type (str): the parent plugin type that is inherited
            from higher level plugins

        """

        # Choosing the correct type according to the plugin name, or its last
        # known authorized parend plugin
        if self.is_allowed(key):
            plg_type = key

        else:
            plg_type = parent_plg_type

        # Initializing the attribute plugin to store type, name and version
        # Update the type if not already given
        plg = getattr(self, "plugin", None)

        if plg is None:
            self.plugin = self.from_dict(
                {"name": None, "version": None, "type": plg_type}
            )

        else:
            if getattr(self.plugin, "type", None) is None:
                self.plugin.type = plg_type

        return plg_type

    @classmethod
    def save_loaded(cls, plg):
        """Saves all loaded plugins to the class

        Args:
            plg (Plugin): plugin to save

        """

        name = plg.plugin.name
        version = plg.plugin.version
        plugin_type = plg.plugin.type

        cls.loaded_instances[cls.plugin_key(name, version, plugin_type)] = plg

    @classmethod
    def _save_refplugins(cls, plg):
        """Save all level 0 attributes to the class.
        Should be called only once

        Args:
            plg (Plugin): plugin to save

        """

        for attr in plg.attributes:
            plg_attr = getattr(plg, attr)

            if issubclass(type(plg_attr), Plugin):
                plg_attr._load_plugin_type(attr)

                name = plg_attr.plugin.name
                version = plg_attr.plugin.version
                plugin_type = plg_attr.plugin.type

                setattr(plg_attr, "isreference", True)

            else:
                name = version = plugin_type = attr

            cls.reference_instances[plugin_type] = plg_attr

    @classmethod
    def _save_subrefplugins(cls, plg, parent_plg_type="Setup", tree=""):
        """Save all attributes to the class.
        Should be called only once

        Args:
            plg (Plugin): plugin to save

        """

        if not hasattr(plg, "attributes"):
            return

        for attr in plg.attributes:
            plg_attr = getattr(plg, attr)

            if issubclass(type(plg_attr), Plugin):
                plg_attr._load_plugin_type(attr, parent_plg_type)

                name = plg_attr.plugin.name
                version = plg_attr.plugin.version
                plugin_type = plg_attr.plugin.type

            else:
                name = version = plugin_type = attr

            plg_tree = "{}/{}".format(tree, attr)
            cls._save_subrefplugins(plg_attr, parent_plg_type, tree=plg_tree)

            if plugin_type in cls.subreference_instances:
                cls.subreference_instances[plugin_type][plg_tree] = plg_attr
            else:
                cls.subreference_instances[plugin_type] = {plg_tree: plg_attr}

    @classmethod
    def from_dict(cls, def_dict, orig_name="", convert_none=False, **kwargs):
        """Loads a recursive dictionary structure into a Plugin

        Args:
            def_dict (dict): the definition dictionary
            orig_dict (dict): the definition dictionary used at level 0

        Returns:
            Plugin
        """
        # Loop over keys to be loaded
        # Recursively initialize from dictionary
        plg = cls()

        for key in def_dict:
            if isinstance(def_dict[key], dict):
                setattr(
                    plg,
                    key,
                    cls.from_dict(
                        def_dict[key],
                        orig_name=key,
                        convert_none=convert_none,
                        **kwargs
                    ),
                )

            # Initializes empty keys as Plugins
            elif def_dict[key] is None and convert_none:
                setattr(plg, key, cls())

            else:
                setattr(plg, key, def_dict[key])

        # Saves the definition keys of the original dictionary into the
        # output plugin
        plg.attributes = list(def_dict.keys())

        # Saves the name of the plugin as specified in the configuration file
        plg.orig_name = orig_name

        return plg

    @classmethod
    def to_dict(cls, plg, exclude_patterns=[], full=False):
        """Turns a Plugin to a dictionary for easier saving.

        Args:
            plg (Plugin): a Plugin instance to be saved as a dictionary
        """

        # If the input has no method 'attributes', just return it
        if not hasattr(plg, "attributes") or isinstance(plg, list):
            return plg

        out = {}

        for attr in dir(plg):
            plg_attr = getattr(plg, attr)

            # Otherwise, save the attribute value
            if (
                not ismethod(plg_attr)
                and not isinstance(plg_attr, FunctionType)
                and not isinstance(plg_attr, ModuleType)
                and not attr[0] == "_"
                and attr not in dir(cls)
                and attr
                not in [
                    "absolute_import",
                    "loaded_class",
                    "loaded_data",
                    "loaded_attributes",
                    "loaded_requirements",
                    "loaded_template",
                    "orig_name",
                    "logfile",
                    "datei",
                    "datef",
                    "workdir",
                    "verbose",
                    "requirements",
                    "attributes",
                    "isreference",
                    "Method_type",
                    "MethodType",
                    "default_requirements",
                    "mapper",
                ]
                + exclude_patterns
                and "Feature" not in str(plg_attr)
            ):

                # If the attribute is a Plugin sub-class,
                # recursively call to_dict
                if hasattr(plg_attr, "to_dict"):
                    out[attr] = cls.to_dict(plg_attr, exclude_patterns, full)

                elif full:
                    if type(plg_attr) in [np.ndarray, pd.core.frame.DataFrame]:
                        out[attr] = "{} {}".format(
                            type(plg_attr), plg_attr.shape
                        )

                    elif type(plg_attr) == dict:
                        out[attr] = "{} {} keys".format(
                            type(plg_attr), len(plg_attr.keys())
                        )

                    elif type(plg_attr) == list:
                        out[attr] = "{} len({})".format(
                            type(plg_attr), len(plg_attr)
                        )
                    else:
                        out[attr] = str(plg_attr)

                else:
                    out[attr] = plg_attr

        # Removing the __class__ attribute
        if "__class__" in list(out.keys()):
            del out["__class__"]

        return out

    @classmethod
    def to_yaml(cls, plg, yaml_file, full=True):
        """Write a Yaml from a loaded plugin"""

        plg_dict = cls.to_dict(plg)

        with open(yaml_file, "w") as f:
            yaml.dump(plg_dict, f)

    @classmethod
    def print_default(cls, plg):
        """Print default values if available"""

        if not hasattr(plg, "default_values"):
            print(
                """{} ({}, {}, {}) has no default values""".format(
                    plg, plg.plugin.name, plg.plugin.version, plg.plugin.type
                )
            )
            return

        print(
            """The default values of {} ({}, {}, {}) are:""".format(
                plg, plg.plugin.name, plg.plugin.version, plg.plugin.type
            )
        )
        for k in plg.default_values:
            print("- {}:\t{}".format(k, plg.default_values[k]))

    @classmethod
    def childclass_factory(cls, plg_orig, child_type=None, overwrite=False):
        """Generates an instance of one of Plugin's child classes. Transfers
        all existing attributes in the argument plugin to the output
        child-class instance

        Args:
            plg_orig (Plugin): the plugin to turn into a child-class instance
            child_type (str): sub-class type to generate if not available in
                              plg_orig.plugin.type
            overwrite (bool): overwrite the class type of the origin plugin

        Return:
            child_plg: a plugin with all the attributes from plg_orig,
            but as a child-class instance
        """

        if getattr(plg_orig.plugin, "type", None) is None:
            if child_type is None:
                raise PluginError(
                    "The Child-class factory was called on a plugin "
                    "that was not correctly initialized: {} / {}".format(
                        plg_orig, plg_orig.orig_name
                    )
                )

            else:
                plg_type = child_type
                plg_orig._load_plugin_type(plg_type)

        else:
            plg_type = plg_orig.plugin.type

        # If the type is not referenced, don't do anything
        if not cls.is_allowed(plg_type):
            return plg_orig

        # Load the subclass
        child_plg = cls.get_subclass(plg_type)(plg_orig=plg_orig)

        return child_plg

    def initiate(self, plg_type=None):
        """Initializes a Plugin, i.e., loads functions from registered
        plugins

        Args:
            plg_type (str): the type of plugin to load; this should
            correspond to one of the defined child-classes

        Return:
            module: a python module as registered in pyCIF

        """

        # It there is no attribute 'plugin', can't initialize anything as
        # python will not know the plugin type, name and version
        plugin = getattr(self, "plugin", None)
        if plugin is None:
            return

        # Load plugin IDs
        name = getattr(plugin, "name", None)
        version = getattr(plugin, "version", None)
        if plg_type is None:
            plg_type = plugin.type

        # Load registered module if the Plugin's name and version are not
        # default empty strings
        if name is not None:
            module = self.get_registered(name, version, plg_type)

            # Attributing all module functions to the plugin
            functions = [
                f
                for f in dir(module)
                if f not in getattr(module, "attributes", [])
            ]
            for attr in functions:
                setattr(self, attr, getattr(module, attr))

            # Ini_data should be define as a Method Type
            if hasattr(module, "ini_data"):
                self.ini_data = MethodType(module.ini_data, self)

            return module
    
    @classmethod
    def flushall(cls):
        cls.subreference_instances = {}
        cls.loaded_instances = {}
        cls.reference_instances = {}