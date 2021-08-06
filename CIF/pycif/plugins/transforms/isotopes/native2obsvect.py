

def native2obsvect(
    transform,
    xmod,
    mapper,
    mod_input,
    ddi,
    ddf,
    mode,
    runsubdir,
    workdir,
    trans_mode,
):

    dfout = xmod.iloc[xmod["indorig"].unique()]
    for isostep in transform.isospecies.attributes[::-1]:
        isostp = getattr(transform.isospecies, isostep)
        names = isostp.inputs.names
        outputs = isostp.outputs.names
        r_std = isostp.inputs.standard
        iso_mass = isostp.outputs.iso_mass
        spec_mass = isostp.outputs.spec_mass
        unit = isostp.unit

        if not mod_input == "obs":
            return xmod

        dfs = [xmod.loc[xmod["parameter"] == s.lower()] for s in outputs]

        mask_iso = dfs[0]["parameter_ref"] == names[0].lower()
        mask_tot = dfs[0]["parameter_ref"] == names[1].lower()

        isotopologue_iso = dfs[1].loc[mask_iso, "sim"]
        spec_ref_iso = dfs[0].loc[mask_iso, "sim"]
        isotopologue_tot = dfs[1].loc[mask_tot, "sim"]
        spec_ref_tot = dfs[0].loc[mask_tot, "sim"]

        # Saving forward simulations for later use by adjoint
        isostp.fwd_data = {
            "isotopologue_iso": isotopologue_iso,
            "spec_ref_iso": spec_ref_iso,
            "isotopologue_tot": isotopologue_tot,
            "spec_ref_tot": spec_ref_tot,
        }

        iso = isotopologue_iso / spec_ref_iso
        if unit == "volume":
            tot = isotopologue_tot + spec_ref_tot

        else:
            tot = (
                spec_ref_tot / iso_mass[0] * spec_mass
                + isotopologue_tot / iso_mass[1] * spec_mass
            )
            iso *= iso_mass[0] / iso_mass[1]

        iso /= r_std
        iso = (iso - 1) * 1000

        # Filling the original datastore with correct values
        # Changing names as well
        ind_iso = dfs[0].loc[mask_iso, "indorig"]
        ind_tot = dfs[0].loc[mask_tot, "indorig"]

        dfout.iloc[ind_iso, dfout.columns.get_loc("sim")] = iso
        dfout.iloc[ind_tot, dfout.columns.get_loc("sim")] = tot

        dfout.iloc[ind_iso, dfout.columns.get_loc("parameter")] = names[
            0
        ].lower()
        dfout.iloc[ind_tot, dfout.columns.get_loc("parameter")] = names[
            1
        ].lower()

        # Applying tangent-linear
        if mode == "tl":
            isotopologue_iso_tl = dfs[1].loc[mask_iso, "sim_tl"]
            spec_ref_iso_tl = dfs[0].loc[mask_iso, "sim_tl"]
            isotopologue_tot_tl = dfs[1].loc[mask_tot, "sim_tl"]
            spec_ref_tot_tl = dfs[0].loc[mask_tot, "sim_tl"]

            iso_tl = (
                isotopologue_iso_tl / spec_ref_iso
                - isotopologue_iso / spec_ref_iso ** 2 * spec_ref_iso_tl
            )
            if unit == "volume":
                tot_tl = isotopologue_tot_tl + spec_ref_tot_tl

            else:
                tot_tl = (
                    isotopologue_iso_tl * iso_mass[0] / spec_mass
                    + spec_ref_tot_tl * iso_mass[1] / spec_mass
                )
                iso_tl *= iso_mass[0] / iso_mass[1]

            iso_tl /= r_std / 1000

            dfout.iloc[ind_iso, dfout.columns.get_loc("sim_tl")] = iso_tl
            dfout.iloc[ind_tot, dfout.columns.get_loc("sim_tl")] = tot_tl

    return dfout
