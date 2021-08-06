Implementation
---------------------

Following Eq. :eq:`bayes` to :eq:`enkf`, operating both variational inversions and EnKFs requires the structure summarized below.
The generic call graph of the CIF could easily be extended to any other data assimilation method, with the observation operator as the common block between all methods.
The observation operator is required to translate information from the control space to the observation space, and reversely in the case of methods using the adjoint.
More specific intermediate mappings are required for interfacing the CIF with the chosen dynamical model.
They include mappings from/to algebraic spaces used by the inversion to/from practical spaces compatible with the dynamical model inputs and outputs.
Such operations depend on the definitions of the control and observation vectors, on the model spatial and temporal resolution,
as well as on the format of input and output files for the model.

Please check :doc:`here <../documentation/plugins/dependencies>` for further details on how the CIF blocks are built and linked with each other.

.. image:: /images/CIFscheme.png
    :target: ../_images/CIFscheme.png