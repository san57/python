from .native2state import native2state
from .state2native import state2native

requirements = {"model": {"any": True, "empty": False}}


def ini_mapper(
    transform, transform_type, inputs={}, outputs={}, backup_comps={}
):
    # Keeps same info from the tracer
    trid = (transform.component[0], transform.parameter[0])

    mapper = {
        "inputs": {trid: dict(inputs.get(trid, {}), **{"force_dump": False})},
        "outputs": {trid: inputs.get(trid, {})},
    }

    return mapper
