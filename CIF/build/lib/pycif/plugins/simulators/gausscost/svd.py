from __future__ import division

import numpy as np
import pandas as pd
from builtins import range


def svd_init(datastore):
    parameters = datastore["parameter"].unique()

    SVD = {}
    for param in parameters:
        ds = datastore.loc[datastore["parameter"] == param]
        stations = ds["station"].unique()
        ds_obs_ref = pd.DataFrame(
            {
                s: ds.loc[ds["station"] == s, "obs"].resample("1D").mean()
                for s in stations
            }
        )

        # Filling NaNs assuming existing patterns can explain missing values
        valid = ~np.isnan(ds_obs_ref)
        mu_hat = np.nanmean(ds_obs_ref, axis=0, keepdims=True)
        ds_obs = np.where(valid, ds_obs_ref, mu_hat)

        for k in range(20):
            U, s, Vh = np.linalg.svd(ds_obs, full_matrices=False)
            ds_obs[~valid] = (U.dot(np.diag(s).dot(Vh)))[~valid]

        SVD[param] = {"U": U, "Vh": Vh, "s": s}

    return SVD


def svd_cost(self, datastore, Jref):
    """Computes the cost function based on the SVD decomposition"""

    j_o = 0
    parameters = datastore["parameter"].unique()
    for param in parameters:
        mask_param = datastore["parameter"] == param

        # Projecting at daily scale
        ds = datastore.loc[mask_param]
        dds = ds.loc[:, "sim"] - ds.loc[:, "obs"]
        stations = ds["station"].unique()
        resamples = {
            s: dds.loc[ds["station"] == s].resample("1D") for s in stations
        }
        ds_delta = pd.DataFrame({s: resamples[s].sum() for s in stations})
        ds_counts = pd.DataFrame({s: resamples[s].size() for s in stations})

        ds_delta /= ds_counts

        mu_hat = np.nanmean(ds_delta, axis=0, keepdims=True)
        valid = ~np.isnan(ds_delta)
        ds_delta = np.where(valid, ds_delta, mu_hat)

        # Fetching SVD vectors
        svd_vect = self.svd_vectors[param]
        U = svd_vect["U"]
        s_obs = svd_vect["s"]

        # SVD and cost function
        Vh_sim = np.dot(U.T, ds_delta)

        # Errors are chosen as inversely proportional to the singular values
        # It is possible to relax constraints on singular vectors from a
        # given index with the argument 'crop_svd' in the Yaml
        errors = 1.0 / s_obs[:, np.newaxis] ** 0.25
        if hasattr(self, "crop_svd"):
            errors[self.crop_svd:] = np.inf

        j_o += np.nansum(Vh_sim ** 2 / errors ** 2)

        # Now get the adjoint for the gradient
        departures = U.dot(Vh_sim / errors ** 2)

        ds.loc[:, "obs_incr"] = np.nan
        for k, stat in enumerate(stations):
            col_stat = ds_counts.columns.get_loc(stat)
            ds_depart = pd.Series(
                departures[:, col_stat], index=ds_counts.index
            )
            mask_stat = ds["station"] == stat
            index = ds.index[mask_stat].floor("D")

            # Propagating mu_hat to adjoint
            nvalid = valid.loc[:, stat].sum()
            ds_depart.loc[valid.loc[:, stat]] += (
                ds_depart.loc[~valid.loc[:, stat]].sum() / nvalid
            )

            # De-aggregating daily scale
            ds.loc[mask_stat, "obs_incr"] = (
                ds_depart.loc[index].values / ds_counts.loc[index, stat].values
            )

        datastore.loc[mask_param, "obs_incr"] = ds.loc[:, "obs_incr"]

    return j_o
