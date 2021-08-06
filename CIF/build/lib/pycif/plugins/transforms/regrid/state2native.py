import os

import numpy as np
import xarray as xr
from scipy.interpolate import RectBivariateSpline

from .reproject import reproject_emissions

try:
    import cPickle as pickle
except ImportError:
    import pickle


def state2native(
    transform,
    data,
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

    xmod = data.datastore

    for trid in mapper["inputs"]:
        wgt_file = mapper["inputs"][trid]["weight_file"]

        # Inputs domain from the present tracer
        nlon_in = xmod[trid]["tracer"].domain.nlon
        zlon_in = xmod[trid]["tracer"].domain.zlon
        zlonc_in = xmod[trid]["tracer"].domain.zlonc
        nlat_in = xmod[trid]["tracer"].domain.nlat
        zlat_in = xmod[trid]["tracer"].domain.zlat
        zlatc_in = xmod[trid]["tracer"].domain.zlatc

        # Outputs domain from the mapper
        is_lbc = mapper["outputs"][trid].get("is_lbc", False)
        domain_out = mapper["outputs"][trid]["domain"]

        xmod[trid]["spec"] = do_regridding(
            transform,
            xmod,
            trid,
            wgt_file,
            nlon_in,
            nlat_in,
            zlon_in,
            zlat_in,
            zlonc_in,
            zlatc_in,
            domain_out,
            is_lbc,
        )

    data.datastore = xmod

    return data


def do_regridding(
    transform,
    xmod,
    trid,
    wgt_file,
    nlon_in,
    nlat_in,
    zlon_in,
    zlat_in,
    zlonc_in,
    zlatc_in,
    domain_out,
    is_lbc,
):
    if is_lbc:
        nlon_out = domain_out.nlon_side
        zlon_out = domain_out.zlon_side
        zlonc_out = domain_out.zlonc_side
        nlat_out = domain_out.nlat_side
        zlat_out = domain_out.zlat_side
        zlatc_out = domain_out.zlatc_side

    else:
        nlon_out = domain_out.nlon
        zlon_out = domain_out.zlon
        zlonc_out = domain_out.zlonc
        nlat_out = domain_out.nlat
        zlat_out = domain_out.zlat
        zlatc_out = domain_out.zlatc

    # Loading weights if available
    loaded_weights = False
    if os.path.isfile(wgt_file):
        with open(wgt_file, "rb") as f:
            weights = pickle.load(f)
            loaded_weights = True

    # Otherwise, creating them
    var_in = xmod[trid]["spec"][0, 0].values
    meshj, meshi = np.meshgrid(range(nlon_in), range(nlat_in))
    if not loaded_weights:
        if transform.method == "mass-conservation":
            weights = reproject_emissions(
                var_in,
                zlonc_in,
                zlatc_in,
                zlonc_out,
                zlatc_out,
                return_weight=True,
            )

        elif transform.method == "bilinear":
            interp = RectBivariateSpline(
                zlon_in[0],
                zlat_in[:, 0],
                meshi.T,
                bbox=[None, None, None, None],
            )
            meshi_out = interp(
                zlon_out.flatten(order="F"),
                zlat_out.flatten(order="F"),
                grid=False,
            )

            interp = RectBivariateSpline(
                zlon_in[0],
                zlat_in[:, 0],
                meshj.T,
                bbox=[None, None, None, None],
            )
            meshj_out = interp(
                zlon_out.flatten(order="F"),
                zlat_out.flatten(order="F"),
                grid=False,
            )

            imin = np.floor(meshi_out).astype(int)
            jmin = np.floor(meshj_out).astype(int)
            alpha_imax = meshi_out - imin
            alpha_jmax = meshj_out - jmin

            weights = [
                (
                    [i, i, i + 1, i + 1],
                    [j, j + 1, j + 1, j],
                    [
                        (1 - ai) * (1 - aj),
                        (1 - ai) * aj,
                        ai * aj,
                        ai * (1 - aj),
                    ],
                )
                for i, j, ai, aj in zip(imin, jmin, alpha_imax, alpha_jmax)
            ]

        else:
            raise Exception(
                "Regrid method not known: {}\n"
                "Please check your Yaml file".format(transform.method)
            )

        with open(wgt_file, "wb") as f:
            pickle.dump(weights, f)
            loaded_weights = True

    # Applying weights
    nlev = len(xmod[trid]["spec"].lev)
    ntimes = len(xmod[trid]["spec"].time)
    var_out = np.empty((ntimes, nlev, nlat_out, nlon_out))
    var_in = xmod[trid]["spec"].values
    for k, wgt in enumerate(weights):
        iout, jout = np.unravel_index(k, (nlat_out, nlon_out), order="F")
        var_out[..., iout, jout] = np.sum(
            var_in[..., wgt[1], wgt[0]] * wgt[2], axis=2
        )

    times = xmod[trid]["spec"].time.values

    return xr.DataArray(
        var_out, coords={"time": times}, dims=("time", "lev", "lat", "lon")
    )
