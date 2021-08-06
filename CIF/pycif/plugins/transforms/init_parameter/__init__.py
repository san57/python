from .native2state import native2state
from .state2native import state2native

requirements = {"model": {"any": True, "empty": False}}


def ini_mapper(
    transform, transform_type, inputs={}, outputs={}, backup_comps={}
):

    if len(transform.component) == 1:
        outputs = [
            (transform.component[0], param) for param in transform.parameter
        ]
    else:
        outputs = [
            (comp, param)
            for (comp, param) in zip(transform.component, transform.parameter)
        ]

    default_dict = {
        "input_dates_parameters": transform.orig_parameter_plg.input_dates,
        "input_files_parameters": transform.orig_parameter_plg.input_files,
        "tracer": transform.orig_parameter_plg,
        "component": transform.orig_component_plg,
    }

    mapper = {
        "inputs": {},
        "outputs": {trid: default_dict for trid in outputs},
    }

    # Forcing reading if need that info later in the transform pipeline
    for trid in mapper["outputs"]:
        if trid in inputs:
            mapper["outputs"][trid]["force_read"] = inputs[trid].get(
                "force_read", False
            )

    return mapper
