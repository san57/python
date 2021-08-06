Yaml configuration
------------------

pyCIF can be entirely set up with a
`Yaml <http://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html>`__
file such as in the example below. The basic idea of Yaml
syntax is that you can easily define a tree structure using ':' and
indentation, which would be automatically interpreted by Python.

Please note that Yaml is flexible with the keys it accepts for booleans.
The following keys will be automatically accepted as True or False:

y / Y / yes / Yes / YES / n / N / no / No / NO / true / True / TRUE /
false / False / FALSE / on / On / ON / off / Off / OFF

Basic concept
^^^^^^^^^^^^^

For instance, the following Yaml syntax:

.. code-block:: yaml

    root:
       branch1: grub
       branch2: [0, 1]
       branch3:
          subbranch0: some text
          subbranch1: 2004-01-01
          subbranch2: T

will be interpreted into the following Python dictionary:

.. code-block:: python

    {
    'root':
       'branch1': 'grub',
       'branch2': [0,1],
       'branch3':
          'subbranch0': 'some text',
          'subbranch1': datetime.datetime(2004, 01, 01),
          'subbranch2': True
    }

which is itself loaded into a pyCIF plugin automatically integrated with
other plugins and dependencies.

Further examples about the Yaml syntax can be found
`here <http://sweetohm.net/article/introduction-yaml.en.html>`__.

pyCIF structure in Yaml
^^^^^^^^^^^^^^^^^^^^^^^

pyCIF is a combination of so-called plugins interacting with each other
through methods and dependencies. Depending on the required set-up it is
possible to change the version of a specifi plugin consistently to your
usage. The main plugins that a regular user will have to adjust are the
following (click on each plugin to see details and compatible set-ups).
All available plugins as implemented in pyCIF are given :doc:`here. <../documentation/plugins/available>`

.. role:: bash(code)
   :language: bash


1.  :bash:`todo_init` (optional):
        a list of modules to be initialized (same
        names as main entries in the Yaml); pyCIF only computes initialization
        in that case, and not the full running mode; it is recommended to first
        generate pyCIF files with this key and check whether they are
        consistent, before running the full pyCIF (by commenting the key)
2.  :doc:`mode <modes/index>`:
        running mode for pyCIF; should
        be one of: ''forward'', ''variational'', ''analytical''
3.  :doc:`model <modes/index>`:
        the numerical model to
        interface in pyCIF; integrated models: ''LMDZ'', ''CHIMERE''
4.  :doc:`measurements <modes/index>`:
        the list of species to integrate in pyCIF and the data provider to know which format
        to parse
5.  :doc:`fluxes <modes/index>`:
        surface fluxes and related species
6.  :doc:`obsvect <modes/index>`:
        definition of the observation vector
7.  :doc:`controlvect <modes/index>`:
        definition of the control vector
8.  :doc:`obsoperator <modes/index>`:
        definition of the observation operator
9.  :doc:`domain <modes/index>`:
        domain definition
10. :doc:`meteo <modes/index>`:
        meteo type and source directory


Please note that plugin dependencies should be consistently fulfilled in
the Yaml file. It implies that if a plugins requires another plugins to
work, the later plugin should be properly defined in the Yaml file. As
detailed on the page of each plugin, some dependencies are set by
default if not set in the Yaml file as an attempt to make the
configuration file as easy-to-read as possible.
For more details on dependencies, please go :doc:`here. <../documentation/plugins/dependencies>`


One can find an example to run forward simulations with a Gaussian toy model below:

.. include:: yaml_example.rst