import numpy as np
import pandas as pd
import xarray as xr

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

    ddi = min(di, df)

    xmod = data.datastore

    for trid in mapper["inputs"]:
        var_in = xmod[trid]["spec"]
        ntimes, nlev, nlat, nlon = var_in.shape

        # Inputs domain from the present tracer
        nlev_in = xmod[trid]["tracer"].domain.nlev
        sigma_a_in = xmod[trid]["tracer"].domain.sigma_a
        sigma_b_in = xmod[trid]["tracer"].domain.sigma_b

        # Outputs domain from the present tracer
        domain_out = mapper["outputs"][trid]["domain"]
        nlev_out = domain_out.nlev
        sigma_a_out = domain_out.sigma_a
        sigma_b_out = domain_out.sigma_b

        if "is_top" in mapper["outputs"][trid]:
            nlev_out = 1
            sigma_a_out = domain_out.sigma_a[-1:]
            sigma_b_out = domain_out.sigma_b[-1:]

        # Use virtual surface pressure to interpolate
        psurf = transform.psurf
        pres_in = sigma_a_in + psurf * sigma_b_in
        pres_out = sigma_a_out + psurf * sigma_b_out

        if "log_interp" in mapper["outputs"][trid]:
            pres_in = np.log10(pres_in)
            pres_out = np.log10(pres_out)

        pres_tmp = np.sort(np.unique(np.append(pres_in, pres_out)))
        df_pres = pd.DataFrame(range(len(pres_in)), index=pres_in)
        df_pres = (
            df_pres.reindex(pres_tmp)
            .interpolate(method="index")
            .reindex(pres_out)
        )

        var_out = np.empty((ntimes, nlev_out, nlat, nlon))
        for k, dd in enumerate(pres_out):
            ind = df_pres.iloc[k, 0]
            dmin = np.floor(ind).astype(int)
            wgt = ind - dmin
            var_out[:, k] = (
                var_in[:, dmin] * (1 - wgt)
                + var_in[:, min(dmin + 1, nlev_in - 1)] * wgt
            )

        xmod[trid]["spec"] = xr.DataArray(
            var_out,
            coords={"time": var_in.time},
            dims=("time", "lev", "lat", "lon"),
        )

    data.datastore = xmod

    return data
