from pycif.utils.classes.transforms import Transform
from .utils import add_default


def init_transform(self, transform_pipe, transform_type=None):

    # If state vector has no transform argument, just initialize as empty
    if not hasattr(transform_pipe, "transform"):
        transform_pipe.transform = Transform.from_dict({})

    # If observations initialize at the root level
    if transform_type == "obs":
        mapper = {}
        for transform in transform_pipe.transform.attributes:
            # Replacing the transform by a transform class
            transf = getattr(transform_pipe.transform, transform)
            if transf is None:
                transf = Transform.from_dict({}, orig_name=transform)
                setattr(transform_pipe.transform, transform, transf)

            # Updating the general mapper, and creates a local one
            transf_mapper = transf.ini_mapper(transform_type)
            transf.mapper = transf_mapper
            mapper[transform] = transf_mapper

        transform_pipe.transform.mapper = mapper
        return

    # Shortcuts on names
    comps = transform_pipe.components.attributes
    components = transform_pipe.components
    transfs = transform_pipe.transform.attributes
    transforms = transform_pipe.transform

    # Adding model transform if not already included
    # And initializes mapper
    mapper = {}
    backup_comps = {}
    all_inputs = {}
    all_outputs = {}
    run_id = "run_model"
    if "run_model" not in transfs:
        yml_dict = {
            "plugin": {
                "name": "run_model",
                "version": "std",
                "type": "transform",
            }
        }
        run_transf, run_id = add_default(transforms, yml_dict, position="last")

    else:
        run_transf = getattr(transforms, "run_model")

    # Updating the general mapper, and creates a local one
    transf_mapper = run_transf.ini_mapper(
        transform_type,
        inputs=all_inputs,
        outputs=all_outputs,
        backup_comps=backup_comps,
    )
    run_transf.mapper = transf_mapper
    mapper[run_id] = transf_mapper

    all_inputs.update(transf_mapper["inputs"])
    all_outputs.update(transf_mapper["outputs"])

    # Loops backwards on available transformations and update inputs/outputs
    # According to model inputs
    itransf = len(transfs)
    while itransf > 0:
        itransf -= 1
        ntransf_loc = 0
        transform = transfs[itransf]

        # Initializes mapper if not already done
        if transform not in mapper:
            # Replacing the transform by a transform class
            # if not already initialized
            transf = getattr(transforms, transform)
            if transf is None:
                transf = Transform.from_dict({}, orig_name=transform)
                setattr(transforms, transform, transf)

            # Updating the general mapper, and creates a local one
            transf_mapper = transf.ini_mapper(
                transform_type,
                inputs=all_inputs,
                outputs=all_outputs,
                backup_comps=backup_comps,
            )
            transf.mapper = transf_mapper
            mapper[transform] = transf_mapper

            all_inputs.update(transf_mapper["inputs"])
            all_outputs.update(transf_mapper["outputs"])

        else:
            transf_mapper = mapper[transform]

        for trid in transf_mapper["inputs"]:
            if transf_mapper["inputs"][trid].get("force_dump", False):
                # Adding a dump transformation before the present transform
                yml_dict = {
                    "plugin": {
                        "name": "dump_parameter",
                        "version": "std",
                        "type": "transform",
                    },
                    "component": [trid[0]],
                    "parameter": [trid[1]],
                }
                add_default(
                    transforms,
                    yml_dict,
                    position="index",
                    index=itransf,
                    mapper=mapper,
                    init=True,
                    inputs=all_inputs,
                    outputs=all_outputs,
                    backup_comps=backup_comps,
                )
                ntransf_loc += 1
                itransf += 1
                continue

    # Add default transformations prior to other transforms
    itransf = 0
    for trid in all_inputs:
        prm = trid[1]
        cmp = trid[0]

        cmp_in = cmp if cmp in comps else backup_comps.get(cmp, None)
        if cmp_in is None:
            continue

        cmp_plg = getattr(components, cmp_in)

        # Fetch parameters
        # If no parameters, handle the component as a whole
        if not hasattr(cmp_plg, "parameters"):
            params = cmp_plg
            parameters = [""]
        else:
            params = cmp_plg.parameters
            parameters = params.attributes[:]

        param = getattr(params, prm, cmp_plg)

        # Temporal re-indexing if any
        if hasattr(param, "time_interpolation"):
            tinterp = param.time_interpolation
            yml_dict = {
                "plugin": {
                    "name": "time_interpolation",
                    "version": "std",
                    "type": "transform",
                },
                "method": getattr(tinterp, "method", "bilinear"),
                "component": [cmp],
                "parameter": [prm],
            }
            add_default(
                transforms,
                yml_dict,
                position="index",
                index=itransf,
                mapper=mapper,
                init=True,
                inputs=all_inputs,
                outputs=all_outputs,
                backup_comps=backup_comps,
            )

        # Vertical re-gridding if any
        if hasattr(param, "vertical_interpolation"):
            vinterp = param.time_interpolation
            yml_dict = {
                "plugin": {
                    "name": "vertical_interpolation",
                    "version": "std",
                    "type": "transform",
                },
                "method": getattr(vinterp, "method", "linear"),
                "psurf": getattr(vinterp, "psurf", 1013),
                "component": [cmp],
                "parameter": [prm],
            }
            add_default(
                transforms,
                yml_dict,
                position="index",
                index=itransf,
                mapper=mapper,
                init=True,
                inputs=all_inputs,
                outputs=all_outputs,
                backup_comps=backup_comps,
            )

        # Regridding if any
        if hasattr(param, "regrid"):
            yml_dict = {
                "plugin": {
                    "name": "regrid",
                    "version": "std",
                    "type": "transform",
                },
                "dir_wgt": getattr(param.regrid, "dir_wgt", ""),
                "method": getattr(param.regrid, "method", "bilinear"),
                "component": [cmp],
                "parameter": [prm],
            }
            add_default(
                transforms,
                yml_dict,
                position="index",
                index=itransf,
                mapper=mapper,
                init=True,
                inputs=all_inputs,
                outputs=all_outputs,
                backup_comps=backup_comps,
            )

        # Rescaling if any
        if hasattr(param, "unit_conversion"):
            yml_dict = {
                "plugin": {
                    "name": "unit_conversion",
                    "version": "std",
                    "type": "transform",
                },
                "scale": getattr(param.unit_conversion, "scale", 1),
                "component": [cmp],
                "parameter": [prm],
            }
            add_default(
                transforms,
                yml_dict,
                position="index",
                index=itransf,
                mapper=mapper,
                init=True,
                inputs=all_inputs,
                outputs=all_outputs,
                backup_comps=backup_comps,
            )

        # Initializing the transformation chain for this parameter
        # Create a random name to identify the current transformation
        yml_dict = {
            "plugin": {
                "name": "init_parameter",
                "version": "std",
                "type": "transform",
                "newplg": True,
            },
            "component": [cmp],
            "parameter": [prm],
            "orig_parameter_plg": param,
            "orig_component_plg": cmp_plg,
        }
        add_default(
            transforms,
            yml_dict,
            position="index",
            index=itransf,
            mapper=mapper,
            init=True,
            inputs=all_inputs,
            outputs=all_outputs,
            backup_comps=backup_comps,
        )

    transform_pipe.transform.mapper = mapper
