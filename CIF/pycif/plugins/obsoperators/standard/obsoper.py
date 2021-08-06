from __future__ import absolute_import

import datetime

from builtins import zip

from pycif.utils import path
from pycif.utils.check import info
from pycif.utils.datastores.dump import dump_datastore, read_datastore
from .check import check_inputs


def obsoper(
    self,
    inputs,
    mode,
    run_id=0,
    datei=datetime.datetime(1979, 1, 1),
    datef=datetime.datetime(2100, 1, 1),
    workdir="./",
    reload_results=False,
    **kwargs
):
    """The observation operator.
    This function maps information from the control space to the observation
    space and conversely depending on the running mode.

    Generates model inputs from the control vector
    inputs(x)

    Turns observations into model compatible extraction points
    i.e. y0 (fwd) or dy* (adj)

    Runs the model to get model extractions
    i.e. Hx (fwd) or H^T.dy* (adj)

    if FWD
    Turns back model-equivalents into observation space
    Hx(native) -> Hx(obs)

    Generates departures and multiply by R-1
    -> Hx - y0, R^(-1).(Hx - y0)

    if ADJ
    Turns native increments to control space increments
    H^T.dy* (native) -> H^T.dy* (control)

    Translates control increments to chi values
    -> B^(1/2) . H^T.dy*

    Args:
        inputs: can be a control vector or an observation vector, depending
                on the mode
        mode (str): the running mode; should be one of 'fwd', 'tl' or 'adj'
        run_id (int): the ID number of the current run; is used to define
                the current sub-directory
        datei, datef (datetime.datetime): start and end dates of the
                simulation window
        workdir (str): the parent directory of the simulation to run
        reload_results (bool): look for results from pre-computed simulation
                if any

    Returns:
        observation or control vector depending on the running mode
    """

    # Check that inputs are consistent with the mode to run
    check_inputs(inputs, mode)

    # Create of sub- working directory for the present run
    rundir = "{}/obsoperator/{}_{:04d}/".format(workdir, mode, run_id)
    path.init_dir(rundir)

    # Create save directory for chaining sub-simulations
    path.init_dir("{}/chain/".format(rundir, run_id))

    # Return results from previous versions if existing
    if reload_results:
        if mode in ["fwd", "tl"]:
            try:
                # Saving the directory for possible later use by the adjoint
                self.model.adj_refdir = rundir

                # Trying reading the monitor if any
                obsvect = self.obsvect
                obsvect.datastore = read_datastore(
                    "{}/monitor.nc".format(rundir)
                )
                return obsvect

            except IOError:
                info(
                    "There is no monitor file to be recovered. "
                    "Compute the full forward simulation"
                )

        elif mode == "adj":
            try:
                statevect = self.statevect
                statevect.load("{}/statevect.pickle".format(rundir))
                return statevect

            except IOError:
                info(
                    "There is no statevect file to be recovered. "
                    "Compute the full adjoint simulation."
                )

    # Initializing modules and variables from the setup
    model = self.model
    statevect = self.statevect
    obsvect = self.obsvect
    subsimu_dates = model.subsimu_dates

    if mode in ["fwd", "tl"]:
        obsvect.datastore.loc[:, "sim"] = 0.0
        obsvect.datastore.loc[:, "sim_tl"] = 0.0

        # Dumps control vector in forward and tl modes
        statevect.dump(
            "{}/statevect.pickle".format(rundir),
            to_netcdf=getattr(statevect, "save_out_netcdf", False),
            dir_netcdf="{}/statevect/".format(rundir),
        )

    elif mode == "adj":
        obsvect = inputs
        statevect.dx = 0 * statevect.x
        subsimu_dates = subsimu_dates[::-1]

    # Loop through model periods and runs the model
    self.missingperiod = False
    for di, df in zip(subsimu_dates[:-1], subsimu_dates[1:]):
        # Create a sub-path for each period
        runsubdir = rundir + min(di, df).strftime("%Y-%m-%d_%H-%M")
        _, created = path.init_dir(runsubdir)

        # If the sub-directory was already created,
        # the observation operator considers that the sub-simulation
        # was already properly run, thus passing to next sub-periods
        # Compute the sub-simulation anyway if some previous periods
        # were missing (as stored in self.missingperiod)
        do_simu = (
            created
            or not getattr(self, "autorestart", False)
            or self.missingperiod
        )
        self.missingperiod = do_simu

        # Some verbose
        info(
            "Running {} for the period {} to {}".format(
                model.plugin.name, di, df
            )
        )
        info("Running mode: {}".format(mode))
        info("Sub-directory: {}".format(runsubdir))

        # Prepare observation vector at model resolution
        model.dataobs = obsvect.obsvect2native(
            di, df, mode, runsubdir, workdir
        )
        model.dataobs = self.do_transforms(
            obsvect.transform,
            model.dataobs,
            obsvect.transform.mapper,
            "obs",
            di,
            df,
            mode,
            runsubdir,
            workdir,
            trans_mode="fwd",
        )

        # If the simulation was already carried out, pass to next steps
        # If a sub-period was missing, following ones will be recomputed even
        # if available
        if do_simu:
            # Writing observations for on-line model extraction if any
            model.native2inputs(model.dataobs, "obs", di, df, runsubdir, mode)

        # Prepare inputs for the model and carry out on the fly transformations
        mapper = statevect.transform.mapper
        self.do_transforms(
            statevect.transform,
            statevect,
            mapper,
            "state",
            di,
            df,
            mode,
            runsubdir,
            workdir,
            trans_mode="fwd",
            do_simu=do_simu,
            onlyinit=getattr(self, "onlyinit", False),
        )

        # If only initializing inputs, continue to next sub-period
        if getattr(self, "onlyinit", False):
            continue

        # Update observation vector if necessary
        if mode in ["fwd", "tl"] and obsvect.datastore.size > 0:
            # Read outputs
            model.outputs2native({}, "obs", di, df, runsubdir, mode)
            # Apply obs transformation and update obs datastore
            model.dataobs = self.do_transforms(
                obsvect.transform,
                model.dataobs,
                obsvect.transform.mapper,
                "obs",
                di,
                df,
                mode,
                runsubdir,
                workdir,
                trans_mode="inv",
            )
            obsvect.native2obsvect(model.dataobs, di, df, runsubdir, workdir)

        # Update state vector if necessary
        elif mode == "adj":
            mapper = statevect.transform.mapper
            self.do_transforms(
                statevect.transform,
                statevect,
                mapper,
                "state",
                di,
                df,
                mode,
                runsubdir,
                workdir,
                trans_mode="inv",
            )

        # Keep in memory the fact that it is (or not) a chained simulation
        model.chain = min(di, df)

    # If only initializing inputs, exit
    if getattr(self, "onlyinit", False):
        info("The run was correctly initialized")
        return

    # Re-initalizing the chain argument
    if hasattr(model, "chain"):
        del model.chain

    # Dump observation vector for later use in fwd and tl modes
    # Otherwise dumps the control vector
    if mode in ["fwd", "tl"]:
        dump_type = obsvect.dump_type
        dump_datastore(
            obsvect.datastore,
            file_monit="{}/monitor.{}".format(rundir, dump_type),
            mode="w",
            dump_type=dump_type,
        )

    elif mode == "adj":
        statevect.dump("{}/statevect.pickle".format(rundir))

    # Cleaning unnecessary files
    if getattr(model, "autoflush", False):
        info("Flushing unnecessary files in {}".format(rundir))
        model.flushrun(rundir, mode)

    # Returning the output object depending on the running mode
    if mode in ["fwd", "tl"]:
        return obsvect

    if mode == "adj":
        return statevect
