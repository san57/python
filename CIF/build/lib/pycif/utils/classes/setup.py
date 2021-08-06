import os

import yaml
from builtins import map

import pycif.utils.dates as dates
from pycif.utils.check import init_log, info
from pycif.utils.check.errclass import PluginError
from pycif.utils.yml import ordered_load
from .baseclass import Plugin


class Setup(Plugin):
    @classmethod
    def run_simu(cls, args):
        # Dealing with relative and variable path
        def_file = os.path.abspath(os.path.expanduser(args["def_file"]))

        # Loading Yaml
        setup = cls.yaml_to_setup(def_file)
        setup = cls.load_config(setup)

        # Copying Yaml file for traceability of simulations
        os.system("cp " + setup.def_file + " " + setup.workdir + "/")

        # Saving the loaded configuration
        if getattr(setup, "dump_config", False):
            cls.to_yaml(
                setup,
                "{}/loaded.{}".format(
                    setup.workdir, os.path.basename(setup.def_file)
                ),
            )

        # Run the mode
        if getattr(getattr(setup, "mode", None), "loaded_requirements", False):
            return setup.mode.execute(**args)

        else:
            info(
                "pycif has correctly been initialized "
                "but no execution mode was specified"
            )

    @classmethod
    def yaml_to_setup(cls, config_file):
        # TODO: Allow for anothere type of files than yaml?
        config_dict = cls.from_yaml(config_file)

        # Looking for others Yaml configs files in the main config_file
        config_dict = cls.yaml_subconfigs(config_dict)

        # Load a dictionary to a Setup recursive object
        setup = cls.from_dict(config_dict, convert_none=True)

        return setup

    @classmethod
    def load_config(cls, setup):

        # Creates and initializes the log file
        logfile, workdir = init_log(
            setup.logfile, setup.workdir, setup.verbose
        )
        cls.config_info(setup)

        setup.logfile = logfile
        setup.workdir = workdir

        # Initialize every plugin, requirements and data
        cls.load_setup(setup, level=0)

        return setup

    @classmethod
    def from_yaml(cls, def_file):
        """Generates a dictionary including all pyCIF parameters

        Args:
            def_file (string) : Path to the definition file
                                Handles both absolute and relative paths

        Returns:
            config_dict (dictionary): Dictionary populated with all pyCIF
            parameters

        """

        yml_file = os.path.abspath(os.path.expanduser(def_file))

        try:
            with open(yml_file, "r") as f:
                config_dict = ordered_load(f)

                config_dict["def_file"] = yml_file

                if "datei" in config_dict:
                    # Converting dates to datetime if necessary
                    config_dict["datei"] = dates.date2datetime(
                        config_dict["datei"]
                    )
                    config_dict["datef"] = dates.date2datetime(
                        config_dict["datef"]
                    )
                return config_dict

        except IOError as e:
            info("Couldn't find config file: {}".format(yml_file))
            info("Please check directories")
            raise e

        except yaml.scanner.ScannerError as e:
            info("Error in the syntax of config file: {}".format(yml_file))
            raise e

    @classmethod
    def yaml_subconfigs(cls, config_dict):
        for key, value in config_dict.items():
            if isinstance(value, dict):
                config_dict[key] = cls.yaml_subconfigs(value)
            else:
                if key == "file_yaml":
                    if not os.path.isfile(value):
                        raise OSError(
                            "The Yaml path given is not a file : "
                            "{}".format(value)
                        )
                    if not os.path.exists(value):
                        raise OSError(
                            "The Yaml path given is not valid "
                            "{}".format(value)
                        )
                    config_dict = cls.from_yaml(value)

        return config_dict

    @classmethod
    def config_info(cls, setup):
        """Prints out main input parameters for pyCIF
        """

        verbose_txt = [
            "pyCIF has been initialized with the following parameters:",
            "Yaml configuration file: {}".format(setup.def_file),
            "Log file: {}".format(setup.logfile),
            "Start date: {}".format(setup.datei),
            "End date: {}".format(setup.datef),
            "Working directory: {}".format(setup.workdir),
        ]

        list(map(lambda v: info(v), verbose_txt))

    @classmethod
    def load_setup(
        cls, plg, parent_plg_type=None, level=999, tree="", **kwargs
    ):
        """Loads a Setup plugin.
        Loops recursively over all attributes of the setup to load:
        1) sub-plugins are initialized as Plugin child-class templates (
        Domain, ObsVect, Model, etc);
        2) instances are saved to the Plugin class to be accessible for
        anywhere later one.

        This allows modifications of the data of a given plugin at some place
        of the code to be automatically forwarded to the rest of the code

        Args:
            self (Setup): the setup to load
            parent_plg_type (str): the last recognized plugin type that is
            inherited by children

        """

        # Update orig_dict if not yet defined
        if level == 0:
            # Saves level 0 entries as reference plugins in requirements
            cls._save_refplugins(plg)
            cls._save_subrefplugins(plg)

        # Initialize default values if any
        if hasattr(plg, "default_values"):
            for k in plg.default_values:
                if not hasattr(plg, k):
                    setattr(plg, k, plg.default_values[k])

        # Loop over self attributes and load them as other Class if necessary
        # If an argument 'todo_init' was specified, initialize only listed plg
        if "todo_init" in cls._get_attribute_list(plg):
            attributes = plg.todo_init

        else:
            attributes = [
                a for a in cls._get_attribute_list(plg) if a != "plugin"
            ]

        # Keep in memory the root plg_type
        root_plg_type = parent_plg_type

        for attr in attributes:
            plg_attr = getattr(plg, attr)
            plg_tree = "{}/{}".format(tree, attr)

            # Re-initializing parent type to the root
            parent_plg_type = root_plg_type

            # For reference instances, check whether the Plugin was already
            # initialized as requirement; if so, just take it from reference
            if (
                attr in cls.reference_instances
                and getattr(plg_attr, "isreference", False)
                and getattr(
                    cls.reference_instances.get(attr, None),
                    "loaded_class",
                    False,
                )
            ):
                setattr(plg, attr, cls.reference_instances[attr])
                continue

            # If not a Plugin, continue
            if not issubclass(type(plg_attr), Plugin):
                continue

            # If is still a Setup class, means that should be processed and
            # Initialized
            if isinstance(plg_attr, Setup) and not getattr(
                plg_attr, "loaded_class", False
            ):
                # Load the plugin type depending on the attribute name
                # Do nothing if the attribute is named 'plugin'
                if attr != "plugin":
                    parent_plg_type = plg_attr._load_plugin_type(
                        attr, parent_plg_type
                    )

                # Build a child sub-class and
                # overwrite the Setup class if needed
                plg_attr = cls.childclass_factory(
                    plg_attr, child_type=parent_plg_type
                )

                # Keep in memory that the current attribute class is loaded
                plg_attr.loaded_class = True

            # Load all attributes recursively if not already done
            if not getattr(plg_attr, "loaded_attributes", False):
                cls.load_setup(
                    plg_attr,
                    parent_plg_type,
                    level=level + 1,
                    tree=plg_tree,
                    **kwargs
                )
                plg_attr.loaded_attributes = True

            # Initializes the plugin from registered module if any
            if hasattr(plg_attr, "initiate_template") and not getattr(
                plg_attr, "loaded_template", False
            ):
                plg_attr.initiate_template()

                # Saves the plugin to the class,
                # so it is accessible by everyone anywhere
                # (including its attributes and stored data)
                if hasattr(plg_attr, "plugin"):
                    name = plg_attr.plugin.name
                    version = plg_attr.plugin.version
                    plg_type = plg_attr.plugin.type

                    if (
                        not cls.is_loaded(name, version, plg_type)
                        and name is not None
                    ):
                        cls.save_loaded(plg_attr)

                plg_attr.loaded_template = True

            # If requirements are not already loaded
            if not getattr(plg_attr, "loaded_requirements", False):
                # Load requirements
                cls._check_requirements(
                    plg_attr, parent_plg_type, level, **kwargs
                )

                # The plugin has been correctly loaded at this point
                plg_attr.loaded_requirements = True

            # Initializes the plugin data
            if hasattr(plg_attr, "ini_data") and not getattr(
                plg_attr, "loaded_data", False
            ):
                plg_attr.ini_data(**kwargs)
                plg_attr.loaded_data = True

            # Linking present plugin to reference level 0 if needed
            if getattr(plg_attr, "isreference", False):
                cls.reference_instances[attr] = plg_attr

            # Updating sub-references if needed
            else:
                if attr not in super(Setup, cls).subreference_instances:
                    super(Setup, cls).subreference_instances[attr] = {
                        plg_tree: plg_attr
                    }
                else:
                    super(Setup, cls).subreference_instances[attr][
                        plg_tree
                    ] = plg_attr

            # Attach plugin to the parent plugin
            setattr(plg, attr, plg_attr)

    @classmethod
    def _check_requirements(
        cls, plg, parent_plg_type=None, level=None, **kwargs
    ):
        """Checking that required modules and plugins are loaded.
        If not, load them.

        Requirements are defined in the __init__.py file of the
        corresponding plugin module.

        Args:
            plg (Plugin): a plugin to initialize

        Notes: Some basic parameters are added as requirements to all plugins;
        These are:
            'datei', 'datef', 'workdir', 'logfile', 'verbose'

        """

        # Dealing with default requirements supposed to be given at level 0
        for key in plg.default_requirements:
            if key not in cls._get_attribute_list(plg):
                if key in cls.reference_instances:
                    setattr(plg, key, cls.reference_instances[key])

                else:
                    raise PluginError(
                        "The default key '{}' is not prescribed"
                        "neither in the plugin {}, nor in the "
                        "level 0 of the configuration file".format(key, plg)
                    )

        # Looping over requirements and including them
        for key in plg.requirements:
            key_req = plg.requirements[key]
            fromany = key_req.get("any", False)
            fromsub = key_req.get("subplug", False)
            preftree = key_req.get("preftree", "")
            empty = key_req.get("empty", False)
            name = key_req.get("name", None)
            version = key_req.get("version", "")
            plg_type = key_req.get("type", key)
            newplg = key_req.get("newplg", False)

            # If not from any plugin, but no default value specified, error
            if not fromany and name is None:
                raise PluginError(
                    "{} needs a specific {}, but none was specified \n"
                    "Please check requirements in your module".format(plg, key)
                )

            # If needs a Plugin explicitly defined,
            # look for it at level 0 of setup, or in children,
            # or in unambiguous level N plugins
            plg_tmp = cls._fetch_requirement(
                plg, key, name, version, plg_type, fromsub, empty, preftree
            )

            # If has a prescribed name
            tmp_plugin = getattr(plg_tmp, "plugin", None)
            tmp_name = getattr(tmp_plugin, "name", None)
            tmp_version = getattr(tmp_plugin, "version", None)
            if (
                getattr(getattr(plg_tmp, "plugin", None), "name", None)
                is not None
                and fromany
            ) or (
                tmp_name is not None
                and name == tmp_name
                and version == tmp_version
                and type(plg_tmp) == cls.get_subclass(plg_type)
                and not fromany
            ):
                plg_out = plg_tmp

            # If a default is defined, load from registered
            elif (name is not None and fromany) or (
                tmp_name is None and not fromany
            ):
                plg_out = cls.load_registered(
                    name, version, plg_type, plg_orig=plg_tmp
                )

            # Otherwise, if accepts empty classes
            elif empty:
                plg_out = plg_tmp

            # Error in the yaml if reaching this point
            else:
                raise PluginError(
                    "{} needs a plugin '{}/{}/{}' and an "
                    "inconsistent one was proposed in the Yaml".format(
                        plg, key, name, version
                    )
                )

            if plg_out is None:
                raise Exception(
                    "{} needs a Plugin '{}' to run properly\n"
                    "there is none in its children nor at the level 0 of "
                    "Yaml\n"
                    "Please check your Yaml".format(plg, key)
                )

            # Keep in memory to initialize a new instance of the plugin or not
            plg.plugin.newplg = newplg
            
            # Adding auxiliary attributes if any
            aux_ids = ["name", "version", "type", "any", "subplug", "preftree",
                       "empty", "newplg"]
            for attr in key_req:
                if attr not in aux_ids:
                    setattr(plg_out, attr, key_req[attr])

            # Attaching the requirement to the parent plugin
            setattr(plg, key, plg_out)

        # Load the requirements if not already done
        cls.load_setup(plg, parent_plg_type, level, **kwargs)

    @classmethod
    def _fetch_requirement(
        cls, plg, key, name, version, plg_type, fromsub, empty, preftree
    ):

        pref_plg = [
            s for s in cls.subreference_instances.get(key, []) if preftree in s
        ]

        # If in children
        if key in cls._get_attribute_list(plg):
            plg_tmp = getattr(plg, key)

        # If not in children but at level 0 of Yaml
        elif key in cls.reference_instances:
            plg_tmp = cls.reference_instances[key]

        # If not at level 0, but no ambiguity
        elif len(cls.subreference_instances.get(key, [])) == 1 and fromsub:
            plg_tmp = list(cls.subreference_instances[key].values())[0]

        # If not at level 0, but ambiguity
        elif (
            len(cls.subreference_instances.get(key, [])) >= 1
            and fromsub
            and len(pref_plg) == 1
        ):
            plg_tmp = cls.subreference_instances[key][pref_plg[0]]

        elif empty:
            if cls.is_registered(name, version, plg_type):
                plg_tmp = cls.load_registered(name, version, plg_type)

            else:
                plg_tmp = cls.get_subclass(plg_type)()

        # Error in the yaml if reaching this point
        else:
            plg_tmp = None
            raise PluginError(
                "{} needs a plugin '{}/{}/{}' and an "
                "inconsistent one was proposed in the Yaml".format(
                    plg, key, name, version
                )
            )

        return plg_tmp
