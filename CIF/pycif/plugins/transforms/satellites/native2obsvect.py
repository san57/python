from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import xarray as xr
from builtins import range
from builtins import zip

from pycif.utils.check import info
from pycif.utils.datastores.dump import dump_datastore, read_datastore
from .apply_AK import apply_ak
from .apply_AK import apply_ak_tl
from .vinterp import vertical_interp


def native2obsvect(
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
    """Aggregate simulations at the grid scale to total columns.
    Re-interpolate the model pressure levels to the satellite averaging kernel
    levels. Average using the averaging kernel formula

    """

    ddi = min(di, df)
    datastore = xmod
    ref_indexes = ~xmod.duplicated(subset=["indorig"])
    y0 = xmod.loc[ref_indexes]

    # Building the extended dataframe
    iq1 = y0["iq1"]
    nblinfo = y0["nblinfo"]
    list_satIDs = iq1.loc[y0["is_satellite"]].unique()

    ds_p = datastore.set_index("indorig")[["pressure", "dp", "airm", "hlay"]]
    for satID in list_satIDs:
        satmask = iq1 == satID
        nobs = np.sum(satmask)
        nblloc = nblinfo.loc[satmask].values - 1
        print("satID", satID)
        # Stacking output datastore into levels * nobs
        native_ind_stack = (
            np.flatnonzero(ref_indexes)[satmask]
            + np.arange(transf.model.domain.nlev)[:, np.newaxis]
        )

        # If all nans in datasim, meaning that the species was not simulated
        sim = datastore.loc[:, "sim"].values[native_ind_stack]
        if not np.any(~np.isnan(sim)):
            continue

        # Grouping all data from this satellite
        datasim = xr.Dataset(
            {
                "pressure": (
                    ["level", "index"],
                    np.log(ds_p["pressure"].values[native_ind_stack]),
                ),
                "dp": (
                    ["level", "index"],
                    ds_p["dp"].values[native_ind_stack],
                ),
                "airm": (
                    ["level", "index"],
                    ds_p["airm"].values[native_ind_stack],
                ),
                "hlay": (
                    ["level", "index"],
                    ds_p["hlay"].values[native_ind_stack],
                ),
                "sim": (["level", "index"], sim),
            },
            coords={
                "index": nblloc,
                "level": np.arange(transf.model.domain.nlev),
            },
        )

        if mode == "tl":
            datasim["sim_tl"] = (
                ["level", "index"],
                datastore.loc[:, "sim_tl"].values[native_ind_stack],
            )

        # convert CHIMERE fields to the correct unit
        # from ppb to molec.cm-2 if the satellite product is a column
        if getattr(transf.available_satellites, satID).product == "column":
            keys = ["sim"] + (["sim_tl"] if mode == "tl" else [])
            for k in keys:
                datasim[k] *= datasim["hlay"] / (1e9 / datasim["airm"])

        # Check whether there is some ak
        file_aks = ddi.strftime(
            "{}/obsvect/satellites/infos_{}%Y%m%d%H%M.nc".format(
                transf.workdir, satID
            )
        )
        print("infos file", file_aks)

        # try:
        if True:
            sat_aks = read_datastore(file_aks).to_xarray()

        # except IOError:
        # Assumes total columns?
        # datastore['qdp'] = datastore['sim'] * datastore['dp']
        # groups = datastore.groupby(['indorig'])
        # y0.loc[:, 'sim'] += \
        #     groups['qdp'].sum().values / groups['dp'].sum().values
        #
        # if 'sim_tl' in datastore:
        #     datastore['qdp'] = datastore['sim_tl'] * datastore['dp']
        #     groups = datastore.groupby(['indorig'])
        #     y0.loc[:, 'sim_tl'] += \
        #         groups['qdp'].sum().values / groups['dp'].sum().values
        #    continue

        aks = sat_aks["ak"][nblloc, -2::-1].T

        if getattr(transf.available_satellites, satID).pressure == "Pa":
            pavgs = sat_aks["pavg"][nblloc, ::-1].T
        else:
            pavgs = 100 * sat_aks["pavg"][nblloc, ::-1].T
        # import code
        # code.interact(local=dict(locals(), **globals()))

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
        qa0lus = sat_aks["qa0lu"][nblloc, -2::-1].T

        # Interpolating simulated values to averaging kernel pressures
        # Doing it by chunk to fastensim_ak = 0. * pavgs_mid the process
        # A single chunk overloads the memory,
        # while too many chunks do not take advantage
        # of scipy automatic parallelisation
        # 50 chunks seems to be fairly efficient
        sim_ak = 0.0 * pavgs_mid
        sim_ak_tl = 0.0 * pavgs_mid

        nchunks = getattr(transf, "nchunks", 50)
        chunks = np.linspace(0, nobs, num=nchunks, dtype=int)
        cropstrato = getattr(
            getattr(transf.available_satellites, satID), "cropstrato", 0
        )
        for k1, k2 in zip(chunks[:-1], chunks[1:]):
            info("Compute chunk for satellite {}: {}-{}".format(satID, k1, k2))

            # Vertical interpolation
            xlow, xhigh, alphalow, alphahigh = vertical_interp(
                datasim["pressure"][:, k1:k2].values,
                pavgs_mid[:, k1:k2].values,
                cropstrato,
            )

            # Applying coefficients
            meshout = np.array(pavgs_mid.shape[0] * [list(range(k2 - k1))])
            sim = datasim["sim"][:, k1:k2].values
            sim_ak[:, k1:k2] = (
                alphalow * sim[xlow, meshout] + alphahigh * sim[xhigh, meshout]
            )

            if mode == "tl":
                sim_tl = datasim["sim_tl"][:, k1:k2].values
                sim_ak_tl[:, k1:k2] = (
                    alphalow * sim_tl[xlow, meshout]
                    + alphahigh * sim_tl[xhigh, meshout]
                )

        # # Correction with the pressure thickness
        # # WARNING: there is an inconsistency in the number of levels
        # sim_ak *= (datasim['dp'] * datasim['sim']).sum(axis=0).values \
        #           / (dpavgs * sim_ak).sum(axis=0).values \
        #           * dpavgs.sum(axis=0).values / datasim['dp'].sum(
        # axis=0).values
        # if 'sim_tl' in datasim:
        #     sim_ak_tl *= \
        #         (datasim['dp'] * datasim['sim_tl']).sum(axis=0).values \
        #         / (dpavgs * sim_ak_tl).sum(axis=0).values \
        #         * dpavgs.sum(axis=0).values / datasim['dp'].sum(axis=0).values

        # Applying aks
        nbformula = getattr(transf.available_satellites, satID).formula
        chosenlevel = getattr(
            getattr(transf.available_satellites, satID), "chosenlev", 0
        )

        print("product:", getattr(transf.available_satellites, satID).product)
        print("nbformula:", nbformula)
        print("chosenlev:", chosenlevel)
        y0.loc[satmask, "sim"] = apply_ak(
            sim_ak, dpavgs, aks.values, nbformula, qa0lus.values, chosenlevel
        )

        if mode == "tl":
            y0.loc[satmask, "sim_tl"] = apply_ak_tl(
                sim_ak_tl,
                dpavgs,
                aks.values,
                nbformula,
                qa0lus.values,
                chosenlevel,
                sim_ak,
            )

    # Save forward datastore for later use by adjoint
    file_monit = ddi.strftime(
        "{}/chain/monit_%Y%m%d%H%M.nc".format(transf.model.adj_refdir)
    )
    dump_datastore(
        datastore,
        file_monit=file_monit,
        dump_default=False,
        col2dump=["pressure", "dp", "indorig", "hlay", "airm"],
        mode="w",
    )

    return y0
