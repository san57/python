from __future__ import absolute_import
from __future__ import division

import numpy as np
import xarray as xr
from builtins import range
from builtins import zip

from pycif.utils.datastores.dump import read_datastore
from .apply_AK import apply_ak_ad
from .vinterp import vertical_interp


def obsvect2native(
    transf,
    xmod,
    mapper,
    mod_input,
    di,
    df,
    mode,
    runsubdir,
    workdir,
    trans_mode,
    **kwargs
):
    """De-aggregate total columns to the model level."""
    ddi = min(di, df)
    y0 = xmod

    # Save where a raw is a satellite for later
    y0.loc[:, "is_satellite"] = y0.loc[:, "level"] < 0

    # Number of levels to extract for satellites
    dlev = np.ones(len(y0), dtype=int)
    dlev[y0["is_satellite"]] = transf.model.domain.nlev

    # Index in the original data of the level-extended dataframe
    native_inds = np.append([0], dlev.cumsum())

    # Output index
    idx = np.zeros((native_inds[-1]), dtype=int)
    idx[native_inds[:-1]] = np.arange(len(y0))
    np.maximum.accumulate(idx, out=idx)
    native_inds = native_inds[:-1]

    # Output dataframe
    datacol = "obs_incr" if mode == "adj" else "obs"
    col2process = [
        "tstep",
        "tstep_glo",
        "i",
        "j",
        "level",
        "dtstep",
        "parameter",
        "duration",
        "is_satellite",
        datacol,
    ]

    df = y0.iloc[idx, y0.columns.get_indexer(col2process)]

    # Levels
    sublevels = np.meshgrid(
        list(range(transf.model.domain.nlev)),
        np.ones(np.where(y0["level"] <= 0)[0].size),
    )[0].flatten()
    df.loc[df["level"] <= 0, "level"] = sublevels

    # Building the extended dataframe
    # iq1 = (np.abs(y0['level']) - np.abs((y0['level'] / 10.)
    #                                    .astype(int) * 10)) \
    #    .astype(int)
    iq1 = y0["station"]
    nblinfo = ((y0["level"].astype(int) - y0["level"]) * 1e7).astype(int)
    list_satIDs = iq1.loc[y0["is_satellite"]].unique()
    print("list_satIDS obsvect2native", list_satIDs)
    # Saving original values for later re-aggregation
    df.loc[:, "indorig"] = idx
    df.loc[:, "iq1"] = iq1[idx]
    df.loc[:, "nblinfo"] = nblinfo[idx]

    # Stop here if not adjoint
    if mode != "adj":
        return df

    # Load pressure coordinates from previous run
    file_monit = ddi.strftime(
        "{}/chain/monit_%Y%m%d%H%M.nc".format(transf.model.adj_refdir)
    )
    fwd_pressure = read_datastore(file_monit).set_index("indorig")

    for satID in list_satIDs:
        satmask = iq1 == satID
        nobs = np.sum(satmask)

        nblloc = nblinfo.loc[satmask].values - 1

        # Getting the vector of increments
        obs_incr = y0.loc[satmask, "obs_incr"]
        sim_ak = y0.loc[satmask, "sim"]
        # If all increments are NaNs, just pass to next satellite
        if not np.any(obs_incr != 0.0):
            continue

        # Get target pressure
        native_ind_stack = (
            native_inds[satmask]
            + np.arange(transf.model.domain.nlev)[:, np.newaxis]
        )
        datasim = xr.Dataset(
            {
                "pressure": (
                    ["level", "index"],
                    np.log(fwd_pressure["pressure"].values[native_ind_stack]),
                ),
                "dp": (
                    ["level", "index"],
                    fwd_pressure["dp"].values[native_ind_stack],
                ),
                "airm": (
                    ["level", "index"],
                    fwd_pressure["airm"].values[native_ind_stack],
                ),
                "hlay": (
                    ["level", "index"],
                    fwd_pressure["hlay"].values[native_ind_stack],
                ),
            },
            coords={
                "index": nblloc,
                "level": np.arange(transf.model.domain.nlev),
            },
        )

        # Getting averaging kernels
        file_aks = ddi.strftime(
            "{}/obsvect/satellites/infos_{}%Y%m%d%H%M.nc".format(
                transf.workdir, satID
            )
        )

        try:
            sat_aks = read_datastore(file_aks).to_xarray()

        except IOError:
            # Assumes total columns?
            # groups = fwd_pressure.groupby(['indorig'])
            # df['obs_incr'] = y0.ix[idx, 'obs_incr'] * fwd_pressure['dp'] \
            #                  / groups['dp'].sum().values[idx]
            continue

        # Defining ak info
        #        aks = sat_aks['ak'][nblloc, ::-1][:,1:].T
        aks = sat_aks["ak"][nblloc, -2::-1].T

        if getattr(transf.available_satellites, satID).pressure == "Pa":
            pavgs = sat_aks["pavg"][nblloc, ::-1].T
        else:
            pavgs = 100 * sat_aks["pavg"][nblloc, ::-1].T

        pavgs = xr.DataArray(
            pavgs,
            coords={"index": nblloc, "level": np.arange(aks.level.size + 1)},
            dims=("level", "index"),
        )
        pavgs_mid = xr.DataArray(
            np.log(0.5 * (pavgs[:-1].values + pavgs[1:].values)),
            coords={"index": nblloc, "level": np.arange(aks.level.size)},
            dims=("level", "index"),
        )
        dpavgs = xr.DataArray(
            np.diff(-pavgs, axis=0),
            coords={"index": nblloc, "level": np.arange(aks.level.size)},
            dims=("level", "index"),
        )

        #        qa0lus = sat_aks['qa0lu'][nblloc, ::-1][:,1:].T
        qa0lus = sat_aks["qa0lu"][nblloc, -2::-1].T

        # Applying aks
        nbformula = getattr(transf.available_satellites, satID).formula
        chosenlevel = getattr(transf.available_satellites, satID).chosenlev
        print("nbformula", nbformula)
        obs_incr = apply_ak_ad(
            sim_ak, dpavgs, aks, nbformula, qa0lus, chosenlevel, obs_incr
        )

        # Adjoint of the log-pressure interpolation
        obs_incr_interp = 0.0 * datasim["pressure"].values

        nchunks = getattr(transf, "nchunks", 50)
        chunks = np.linspace(0, nobs, num=nchunks, dtype=int)
        cropstrato = getattr(
            getattr(transf.available_satellites, satID), "cropstrato", False
        )
        for k1, k2 in zip(chunks[:-1], chunks[1:]):
            # Vertical interpolation
            xlow, xhigh, alphalow, alphahigh = vertical_interp(
                datasim["pressure"][:, k1:k2].values,
                pavgs_mid[:, k1:k2].values,
                cropstrato,
            )

            # Applying coefficients
            # WARNING: There might be repeated indexes in a given column
            # To deal with repeated index, np.add.at is recommended
            levmeshout = np.array(
                (k2 - k1) * [list(range(pavgs_mid.shape[0]))]
            ).T
            meshout = np.array(pavgs_mid.shape[0] * [list(range(k1, k2))])

            np.add.at(
                obs_incr_interp,
                (xlow, meshout),
                obs_incr[levmeshout, meshout] * alphalow,
            )

            np.add.at(
                obs_incr_interp,
                (xhigh, meshout),
                obs_incr[levmeshout, meshout] * alphahigh,
            )

        # # Correction with the pressure thickness
        # obs_incr_interp *=\
        #     (dpres * obs_incr_interp).sum(axis=0) \
        #     / (dpavgs * obs_incr).sum(axis=0) \
        #     * dpavgs.sum(axis=0) / dpres.sum(axis=0)

        # convert CHIMERE fields to the correct unit
        # from ppb to molec.cm-2 if the satellite product is a column
        if getattr(transf.available_satellites, satID).product == "column":
            obs_incr_interp *= datasim["hlay"].values / (
                1e9 / datasim["airm"].values
            )

        # Applying increments to the flattened datastore
        df.iloc[
            native_ind_stack.flatten(), df.columns.get_loc("obs_incr")
        ] = obs_incr_interp.flatten()

    return df
