###################################
Standard pyCIF control vector
###################################



.. role:: bash(code)
   :language: bash

The control vector plugin includes sub-routines to compute operation
related to uncertainty matrices, and stores the control vector meta-data
and data themselves. Only a standard control vector including most
commonly used control vector shapes is implemented yet.

Configuration
-------------

The control vector is defined with the following parameters:

- :bash:`components`:
    the different types of components in the control
    vector; accepts: :bash:`fluxes`, :bash:`prescrconcs`, :bash:`prodloss3d` and :bash:`inicond`

    - :bash:`fluxes`: fluxes to be optimized
    - :bash:`species`: each species is defined through a paragraph with the following parameters:

        - :bash:`hresol`: the horizontal resolution; should be one of:
        - :bash:`hpixels`: for individual pixels
        - :bash:`bands`: for zonal and meridional bands
        - :bash:`regions`: for pre-defined regions
        - :bash:`type`: (optional) the type of increments to deal with; should be one of:
            - :bash:`scalar`: (default) increments are applied to the scaling factor of the prior
            - :bash:`physical`: valid only with :bash:`component`=`fluxes` and :bash:`hresol`=`hpixels`
        - :bash:`errtype`: (optional) the type of error; if not specified, a scalar is applied to the prior value; if :bash:`max` is given, computes the max of the neighboring cells (spatially and temporally)
        - :bash:`err`: error as a proportion of prior fluxes
        - :bash:`dir`: directory where data files are stored
        - :bash:`file`: file format of the data
        - :bash:`hcorrelations`: (optional) information on horizontal
            correlations if any
        - :bash:`filelsm`: path to the land-sea mask file
        - :bash:`dump_hcorr`: (optional) dump correlation matrix if True;
            default is True
        - :bash:`dircorrel`: directory where to save the correlation matrix
        - :bash:`sigma_land`: correlation distance over land in km
        - :bash:`sigma_sea`: correlation distance over sea in km
        - :bash:`tcorrelations` : (optional) information on temporal
            correlations if any
        - :bash:`sigma_t`: correlation period in hours
        - :bash:`dump_tcorr`: (optional) dump correlation matrix if True;
            default is True
        - :bash:`dircorrel`: directory where to save the correlation matrix
