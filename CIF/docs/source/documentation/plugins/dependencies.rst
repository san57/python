########################################
Requirements and dependencies in pyCIF
########################################

.. role:: bash(code)
   :language: bash

pyCIF automatically links plugins with each other depending on requirements specified
in the source of the plugins and on elements of the configuration file.

For instance, the observation operator needs to run the numerical model at some point.
For this reason, the :bash:`model` plugin that will be run is attached to the
:bash:`obsoper` plugin.
In the code of :bash:`obsoper`, the model can thus simply be called as follows:

.. code-block:: python

    def obsoper(self):

        # Running the model
        self.model.run(**args, **kwargs)

Another example is a :bash:`model` plugin that needs information about its corresponding :bash:`domain` plugin:

.. code-block:: python

    def some-model-method(self):

        # Fetching domain information
        nlon = self.domain.nlon
        nlat = self.domain.nlat


Types of dependencies
-----------------------

There are two main tpes of dependencies in pyCIF:

1.  :bash:`plugin A` calls methods from :bash:`plugin B`
        For that reason, methods must have a standardized format for arguments and outputs,
        as specified in corresponding pages of the documentation.
2.  :bash:`plugin A` needs data from :bash:`plugin B`
        Similarly, data format needs to follow a specified standardized format

Below is an example graph of dependencies corresponding to the pyCIF configuration file :doc:`here</configuration/yaml_example>`:

.. graphviz:: dependencies.dot
    :align: center


In this example, blue boxes are :bash:`plugins` that are explicitly defined in the configuration file.
Red boxes are implicitly deduced from default values as specified in individual :bash:`plugins`
(see :doc:`here<available>`).
Black arrows stand for method dependencies, while red ones are data dependencies.

Thus, for instance, the :bash:`obsoper` plugin is required by the :bash:`mode` plugin while it is not specified in the Yaml file.
pyCIF automatically initializes the default dependency: :bash:`obsoper` (standard, std).

Behaviour with dependencies
---------------------------

A given plugin might need to use another plugin in different ways, which will determine how the Yaml configuration
should be written and how missing requirements will be dealt with.
Here are the possible behaviours accepted in pyCIF:

* any plugin of correct type fits:




