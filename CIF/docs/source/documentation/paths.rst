############
Paths in Yaml
############

.. role:: bash(code)
   :language: bash

When setting up pyCIF, please prefer absolute paths as much as possible.
However, for conveniency, the YAML parser in pyCIF has been designed to understand environment variables, as well as the tilde expansion.

Below is a simple example of how it works:

.. code-block:: yaml
    :linenos:

    myrefdir: ~/somedir/somesubdir

    myrelativedirfromenvvariable: ${ENV_VARIABLE}/somesubdir

Please note that environment variables must be necessary between :bash:`{}`.

