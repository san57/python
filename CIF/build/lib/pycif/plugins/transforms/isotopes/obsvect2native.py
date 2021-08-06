import numpy as np
import pandas as pd


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
):

    xmod.loc[:, "indorig"] = np.arange(len(xmod))
    xmod.loc[:, "parameter_ref"] = xmod["parameter"]

    isospecies = transf.isospecies.attributes
    if mode == "adj":
        isospecies = isospecies[::-1]

    for isostep in isospecies:
        isostp = getattr(transf.isospecies, isostep)
        names = isostp.inputs.names
        outputs = isostp.outputs.names
        r_std = isostp.inputs.standard
        iso_mass = isostp.outputs.iso_mass
        spec_mass = isostp.outputs.spec_mass
        unit = isostp.unit

        # Just splitting into two sub-species in the dataframe
        mask = xmod["parameter"].str.lower().isin([s.lower() for s in names])
        df = xmod.loc[mask].copy()
        df = pd.concat([df.assign(parameter=iso.lower()) for iso in outputs])

        # Computing adjoint if needed
        if mode == "adj":
            mask_tot = df["parameter_ref"] == names[1].lower()
            mask_iso = df["parameter_ref"] == names[0].lower()
            mask_spec_ref = (df["parameter"] == outputs[0].lower()) & mask_iso
            mask_isotopologue = (
                df["parameter"] == outputs[1].lower()
            ) & mask_iso

            incr = 0 * df["obs_incr"]
            incr_spec_ref = df.loc[mask_spec_ref, "obs_incr"]
            incr_isotopologue = df.loc[mask_isotopologue, "obs_incr"]

            # Fetch simulation from previous forward
            isotopologue = isostp.fwd_data["isotopologue_iso"]
            spec_ref = isostp.fwd_data["spec_ref_iso"]

            incr.loc[mask_tot] = df.loc[mask_tot, "obs_incr"]
            incr.loc[mask_spec_ref] = (
                -isotopologue / spec_ref ** 2 * incr_spec_ref * 1000 / r_std
            )
            incr.loc[mask_isotopologue] = (
                1 / spec_ref * incr_isotopologue * 1000 / r_std
            )

            df.loc[:, "obs_incr"] = incr

        xmod = pd.concat([df, xmod.loc[~mask]])

    return xmod
