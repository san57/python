import copy


def state2native(
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
    **kwargs
):

    xmod = data.datastore

    for isostep in transform.isospecies.attributes:
        isostp = getattr(transform.isospecies, isostep)
        names = isostp.inputs.names
        outputs = isostp.outputs.names
        r_std = isostp.inputs.standard
        iso_mass = isostp.outputs.iso_mass
        spec_mass = isostp.outputs.spec_mass
        unit = isostp.unit

        in_index = []

        in_types = mapper["input_component"]
        out_types = mapper["output_component"]

        keys = ["spec"]
        if mode == "tl":
            keys.append("incr")

        xmod[(outputs[0]), out_types[0]] = copy.deepcopy(
            xmod[(names[0], in_types[0])]
        )
        xmod[(outputs[1]), out_types[1]] = copy.deepcopy(
            xmod[(names[0], in_types[0])]
        )

        signature = xmod[(names[0], in_types[0])]["spec"]
        spec_data = xmod[(names[1], in_types[1])]["spec"]

        # Save signature and data for later use by adjoint
        isostp.fwd_data = {
            "signature": signature.values,
            "spec_data": spec_data.values,
        }

        a_factor = (1 + signature / 1000) * r_std

        xmod[(outputs[0]), out_types[0]]["spec"] = spec_data / (1 + a_factor)
        xmod[(outputs[1]), out_types[1]]["spec"] = (
            spec_data * a_factor / (1 + a_factor)
        )

        # Applying tangent linear
        if mode == "tl":
            spec_data_tl = xmod[(names[1], in_types[1])]["incr"]
            signature_tl = xmod[(names[0], in_types[0])]["incr"]

            xmod[(outputs[0]), out_types[0]]["incr"] = spec_data_tl / (
                1 + a_factor
            ) - spec_data * r_std * signature_tl / (1000 * (1 + a_factor) ** 2)
            xmod[(outputs[1]), out_types[1]][
                "incr"
            ] = spec_data_tl * a_factor / (
                1 + a_factor
            ) + spec_data * r_std * signature_tl / (
                1000 * (1 + a_factor) ** 2
            )

        # Mass correction if units are in mass
        if unit == "mass":
            xmod[(outputs[0]), out_types[0]]["spec"] *= iso_mass[0] / spec_mass
            xmod[(outputs[1]), out_types[1]]["spec"] *= iso_mass[1] / spec_mass

            if mode == "tl":
                xmod[(outputs[0]), out_types[0]]["incr"] *= (
                    iso_mass[0] / spec_mass
                )
                xmod[(outputs[1]), out_types[1]]["incr"] *= (
                    iso_mass[1] / spec_mass
                )

        del xmod[(names[0], in_types[0])]
        del xmod[(names[1], in_types[1])]

    data.datastore = xmod
    return data
