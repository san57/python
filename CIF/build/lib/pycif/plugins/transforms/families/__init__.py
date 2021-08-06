from .native2state import native2state
from .state2native import state2native


def ini_mapper(transform, transform_type, inputs={}, **kwargs):

    parameters_in = transform.parameters_in
    parameter_out = transform.parameter_out
    component = transform.component

    trid_out = (component, parameter_out)

    loc_outputs = {
        trid_out: dict(inputs.get(trid_out, {}), **{"force_dump": False})
    }
    loc_inputs = {
        (component, param): dict(
            inputs.get(trid_out, {}),
            **{"force_dump": False, "force_read": True}
        )
        for param in parameters_in
    }

    mapper = {"inputs": loc_inputs, "outputs": loc_outputs}

    if trid_out in inputs:
        del inputs[trid_out]

    return mapper
