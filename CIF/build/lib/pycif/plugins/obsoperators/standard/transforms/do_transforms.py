from pycif.utils.check import debug


def do_transforms(
    transforms,
    data,
    mapper,
    mod_input,
    di,
    df,
    mode,
    runsubdir,
    workdir,
    trans_mode="fwd",
    do_simu=True,
    onlyinit=False,
    **kwargs
):

    # Reversing dates if adjoint run
    ddi = min(di, df)
    ddf = max(di, df)

    # Type of transform
    isobs = mod_input == "obs"

    if isobs:
        # If empty datastore, do nothing
        if data.size == 0:
            return data

    # Re-ordering the order of tranformations if not in fwd mode
    transf_attributes = transforms.attributes[:]
    if trans_mode != "fwd":
        transf_attributes = transf_attributes[::-1]

    for transform in transf_attributes:
        transf_mapper = mapper[transform]
        transf = getattr(transforms, transform)

        if isobs:
            if trans_mode == "fwd":
                data = transf.obsvect2native(
                    transf,
                    data,
                    mapper,
                    mod_input,
                    ddi,
                    ddf,
                    mode,
                    runsubdir,
                    workdir,
                    trans_mode,
                    do_simu=do_simu,
                    onlyinit=onlyinit,
                )
            else:
                data = transf.native2obsvect(
                    transf,
                    data,
                    mapper,
                    mod_input,
                    ddi,
                    ddf,
                    mode,
                    runsubdir,
                    workdir,
                    trans_mode,
                    do_simu=do_simu,
                    onlyinit=onlyinit,
                )

        else:
            debug(
                [
                    "Doing transform {}: {}".format(
                        transform, transf.plugin.name
                    ),
                    "From inputs: {}".format(transf_mapper["inputs"].keys()),
                    "To outputs: {}".format(transf_mapper["outputs"].keys()),
                ]
            )

            if trans_mode == "fwd":
                data = transf.state2native(
                    transf,
                    data,
                    mapper[transform],
                    mod_input,
                    ddi,
                    ddf,
                    mode,
                    runsubdir,
                    workdir,
                    trans_mode,
                    do_simu=do_simu,
                    onlyinit=onlyinit,
                )
            else:
                data = transf.native2state(
                    transf,
                    data,
                    mapper[transform],
                    mod_input,
                    ddi,
                    ddf,
                    mode,
                    runsubdir,
                    workdir,
                    trans_mode,
                    do_simu=do_simu,
                    onlyinit=onlyinit,
                )

    return data
