import copy


def native2state(
    transform,
    data,
    mapper,
    mod_input,
    ddi,
    ddf,
    mode,
    runsubdir,
    workdir,
    trans_mode,
):

    if mode != "adj":
        return data

    xmod = data["adj_out"]

    for isostep in transform.isospecies.attributes[::-1]:
        isostp = getattr(transform.isospecies, isostep)
        names = isostp.inputs.names
        outputs = isostp.outputs.names
        r_std = isostp.inputs.standard
        iso_mass = isostp.outputs.iso_mass
        spec_mass = isostp.outputs.spec_mass
        unit = isostp.unit

        in_types = isostp.mapper["input_component"]
        out_types = isostp.mapper["output_component"]

        xmod[(names[0], in_types[0])] = copy.deepcopy(
            xmod[(outputs[0], out_types[0])]
        )
        xmod[(names[1], in_types[1])] = copy.deepcopy(
            xmod[(outputs[1], out_types[1])]
        )

        signature = isostp.fwd_data["signature"]
        spec_data = isostp.fwd_data["spec_data"]
        isotopologue_adj = xmod[(outputs[1], out_types[1])]["data"]
        spec_ref_adj = xmod[(outputs[0], out_types[0])]["data"]

        a_factor = (1 + signature / 1000) * r_std

        # Mass correction if units are in mass
        if unit == "mass":
            isotopologue_adj *= iso_mass[1] / spec_mass
            spec_ref_adj *= iso_mass[0] / spec_mass

        # Sensitivity to total
        xmod[(names[1], in_types[1])]["data"] = (
            a_factor * isotopologue_adj + spec_ref_adj
        ) / (1 + a_factor)

        # Sensitivity to signature
        xmod[(names[0], in_types[0])]["data"] = (
            r_std
            * spec_data
            * (isotopologue_adj - spec_ref_adj)
            / 1000
            / (1 + a_factor) ** 2
        )

        del xmod[(outputs[0], out_types[0])]
        del xmod[(outputs[1], out_types[1])]

    data["adj_out"] = xmod

    return data
