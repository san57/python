LMDZ-Dispersion-SACS
--------------------


.. role:: bash(code)
   :language: bash

LMDZ-Dispersion is the offline version of the Global Circulation Model
LMDZ. It was originally customized to be interfaced with `PYVAR <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2005JD006390>`__.
The PYVAR interface was adapted to be implemented to the CIF.
More details can be found :doc:`here</documentation/plugins/models/lmdz>`.

pyCIF requires the following arguments to set a LMDZ simulation up:

-  :bash:`periods`:
    length of simulation sub-segments. Should be compatible
    with Pandas frequency syntax (details
    `here <https://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__);
    this argument should be compatible with LMDZ inputs (most of the time
    monthly sub-segments, which corresponds to :bash:`periods` = '1MS',
    meaning every 1st of the month in pandas); optional; default is '1MS'

-  :bash:`fileexec`:
    the model executable, as compiled in the pyCIF directory

-  :bash:`filedef`:
    a template for the definition file; includes some
    information about the model configuration

-  :bash:`restart`:
    initial concentrations to use for the run:

       -  :bash:`dir`: directory where to find the initial conditions
       -  :bash:`file`: file name format, following :doc:`pyCIF guidelines <../file_formats>`

-  :bash:`chemistry`:
    arguments to handle chemistry with the SACS chemistry
    module (`Pison et al., 2009 <https://doi.org/10.5194/acp-9-5281-2009>`__); optional

        -  :bash:`kinetic`:
            files to fetch variables related to chemical kinetics
            (temperature, pressure, reaction rates)
        -  :bash:`deposition`:
            species with deposition velocity and corresponding
            :doc:`files <../file_formats>`
        -  :bash:`prodloss2d`:
            species with 3D chemical production or destruction
            and corresponding :doc:`files <../file_formats>`
        -  :bash:`prescribed`:
            species with fixed fields; these species are not
            transported; LMDZ need fields as
            molecules.cm :sup:`-3`; it is possible to prescribe
            OH in VMR while setting :bash:`convOH` to :bash:`True` to compute the
            conversion on-the-fly

-  :bash:`species`:
    species to be transported and be applied chemistry in
    the model; dir/file are optional if fluxes for the corresponding
    species are already specified in the control vector, but they are
    necessary if no flux component correspond to that species in the
    control vector

        -  :bash:`restart_id`:
            corresponding IDs in the restart files; these IDs
            are supposed to be in the form :bash:`qXXX`
        -  :bash:`dir`:
            directory to the flux files corresponding to the species
        -  :bash:`file`:
            file format for flux files

-  :bash:`dump`:
    :bash:`True` to save concentration output fields as NetCDF
    files; optional; default is :bash:`False`

Extra arguments below are for expert users only:

-  :bash:`physic`/:bash:`thermals`:
    include physics and thermals.

Below an exemple of a corresponding Yaml paragraph:

.. code-block:: yaml

    model :
      plugin:
        name    : LMDZ
        version : std

      # Length of simulation sub-periods (use a Pandas frequency syntax)
      periods: 1MS

      # Executable
      fileexec : ~/CIF/model_sources/DISPERSION_gch/dispersion.e

      # Definition file (includes some parameters for the simulation)
      filedef  : ~/CIF/model_sources/DISPERSION_gch/def/run.def

      # Initial conditions recovered from:
      restart:
        dir : ~/RESTART/LMDZ/39L/
        file : lmdz5.inca.restart.an%Y.m%mj01.nc

      # Include physics and thermals
      physic : True
      thermals : False

      # Chemistry
      chemistry :
        kinetic :
          dir : ~/In/lmdzinca/
          file : TransCom.vmr.1_scaler.CH2O.m%m.nc
        deposition :
          CH2O :
            dir : ~/In/lmdzinca/
            file : TransCom.vmr.1_scaler.CH2O.m%m.nc
        prodloss3d :
          CH2O :
            dir : ~/In/lmdzinca/
            file : TransCom.vmr.1_scaler.CH2O.m%m.nc
        prescribed :
          OH :
            dir : ~/In/lmdzinca/
            file : TransCom.vmr.1_scaler.CH2O.m%m.nc
            convOH : False
          O1D :
            dir : ~/In/lmdzinca/
            file : TransCom.vmr.1_scaler.CH2O.m%m.nc

      # Species to be transported in the model
        species:
          CH4:
            restart_id: 27
            dir: ~/In/priorflux/flx_monthly/zb/TOTAL/
            file: sflx_CO_CH4_MCF_lmdz9696_%Y_phy.nc
          MCF:
            restart_id: 06
            dir: ~/In/priorflux/flx_monthly/zb/TOTAL/
            file: sflx_CO_CH4_MCF_lmdz9696_%Y_phy.nc
          CO:
            restart_id: 28
            dir: ~/In/priorflux/flx_monthly/zb/TOTAL/
            file: sflx_CO_CH4_MCF_lmdz9696_%Y_phy.nc
          CH2O:
            restart_id: 33
            dir: ~/In/priorflux/flx_monthly/zb/TOTAL/
            file: sflx_CO_CH4_MCF_lmdz9696_%Y_phy.nc

      # Dump outputs into a NetCDF file
      dump: True
