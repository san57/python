################################################
How to add, register and initialize a new plugin
################################################

.. role:: bash(code)
   :language: bash

pyCIF is organized around so-called :bash:`plugins` interacting with each other.
Further details about :bash:`plugins` are given  :doc:`here<../../documentation/plugins/index>`.
Here, you will learn how to add a new :bash:`plugins` of any type.

Adding a :bash:`plugin` to pyCIF
---------------------------------

A :bash:`plugin` is simply a python module.
To add a new plugin of a given class to pyCIF, you need to create a new folder in the correct class folder in
:bash:`pycif/plugins/`.
You then need to create an empty file called :bash:`__init__.py`, so python interprets the new folder as a python module.

.. code-block:: bash

    cd pycif/plugins
    cd the-class-you-need-to-create
    mkdir your-new-plugin
    cd your-new-plugin
    touch __init__.py


Registering your new :bash:`plugin` to pyCIF
--------------------------------------------

:bash:`plugins` in pyCIF are identified with:

    - a type
    - a name
    - a version

When the new plugin is created, it must be registered, so it can be called by other routines of pyCIF.
This can be done through the :bash:`register.py` script located at the class level of :bash:`pycif/plugins/`.

Once the empty plugin is created, edit the :bash:`register.py` script by adding the following lines:

.. code-block:: python

    import your-new-model
    Model.register_plugin('Your-Model-Name', 'Some-version', your-new-model)

The example above is given for a :bash:`model` plugin.
You need to slightly adapt the second line of the example according to the class you are creating.

You can check that your :bash:`plugin` is correctly registered by using the following commands in python:

.. code-block:: python

    from pycif.utils.classes.baseclass import Plugin
    Plugin.print_registered()

Adding requirements to your :bash:`plugin`
------------------------------------------

All :bash:`plugins` are not stand-alone modules, but instead need attributes, data or functions from other plugins.
pyCIF allows you to easily interface :bash:`plugins` with each other through so-called :bash:`requirements`.
Further details on requirements can be found :doc:`here<../../documentation/plugins/dependencies>`.

Information about :bash:`requirements` are specified in the :bash:`__init__.py` file of your :bash:`plugin`.
Below is an example of requirements for the :bash:`model` CHIMERE:

.. code-block:: python

    requirements = {'domain': {'name': 'CHIMERE', 'version': 'std',
                               'empty': False, 'any': False},
                    'chemistry': {'name': 'CHIMERE', 'version': 'gasJtab',
                                  'empty': False, 'any': False},
                    'fluxes': {'name': 'CHIMERE', 'version': 'AEMISSIONS',
                               'empty': True, 'any': False},
                    'meteo': {'name': 'CHIMERE', 'version': 'std',
                               'empty': False, 'any': False}}

With the example above, the :bash:`model` plugin, CHIMERE, requires a :bash:`domain` plugin to work.
In the sub-routines of the :bash:`model` plugin, data from the :bash:`domain` can simply be summoned with:

.. code-block:: python

    def some_model_function(self, *args, **kwargs):

        # Needs the number of grid cells
        nlon = self.domain.nlon
        nlat = self.domain.nlat

The content of the :bash:`requirements` python dictionary is interpreted when a pyCIF run is initializing.
Requirements are fulfilled with regards to the Yaml configuration file.
Depending the values of the arguments, :bash:`name`, :bash:`version`, :bash:`empty` :bash:`any`, pyCIF will need the
:bash:`requirement` to be explicitly specified or not, default :bash:`plugins` could be used, etc.
All details on how this works are given :doc:`here<../../documentation/plugins/dependencies>`.

Initializing your :bash:`plugin`
--------------------------------

Your :bash:`plugin` may need some special operations at initialization.
For instance, it may need to copy files, to compile scripts, to read data, to specified default values, etc.
All these operations can be specified in the :bash:`__init__.py` by creating a function called :bash:`ini_data`.

Below is the :bash:`ini_data` method for the :bash:`model` plugin, CHIMERE:

.. code-block:: python

    def ini_data(plugin, **kwargs):
        """Initializes CHIMERE

        Args:
            plugin (dict): dictionary defining the plugin
            **kwargs (dictionary): possible extra parameters

        Returns:
            loaded plugin and directory with executable

        """

        info("Initializing the model")

        workdir = getattr(plugin, 'workdir', './')

        # Initializes the directory
        path.init_dir('{}/model'.format(workdir))

        # Copying the executables
        target = '{}/model/'.format(workdir)

        source = '{}/src/fwdchimere.e'.format(plugin.direxec)
        shutil.copy(source, target)

        source = '{}/src_tl/tlchimere.e'.format(plugin.direxec)
        shutil.copy(source, target)

        source = '{}/src_ad/achimere.e'.format(plugin.direxec)
        shutil.copy(source, target)

        # Required inputs for running a CHIMERE simulation
        plugin.required_inputs = ['exe', 'nml', 'fluxes',
                                  'meteo', 'inicond', 'boundcond']

        # Default values:
        # period: '1D'
        plugin.periods = getattr(plugin, 'periods', '1D')

        # Number of hours per period
        plugin.nhours = int(pd.to_timedelta(plugin.periods).total_seconds() / 3600)
        plugin.nho = '{:.0f}'.format(plugin.nhours)

        # Replace name for AEMISSION files
        plugin.fluxes.file = plugin.fluxes.file.format(nho=plugin.nho)

        return plugin

During initialization of the :bash:`model` plugin, CHIMERE, the following operations are carried out:

    * creating a new directory :bash:`model/` in the working directory with CHIMERE executables to be called later for simulations
    * specifying a default value for the attribute :bash:`period` if not given in the configuration file
    * modifying the template name of CHIMERE fluxes


Attaching functions to the plugin
---------------------------------

Any plugin is a combination of functions that will be called by itself or by other plugins.
It is necessary to attach functions to your plugin so that they can be called this way in other scripts:

.. code-block:: python

    my-plugin.my-function(*args, **kwargs)

To do so, you simply need to import the function of interest at the module level, i.e., in the file :bash:`__init__.py`.
For any function of interest the import line should be added to the :bash:`__init__.py` file.

.. code-block:: python

    import my-function














