from .native2state import native2state
from .state2native import state2native

requirements = {"model": {"any": True, "empty": False}}


def ini_mapper(
    transform, transform_type, inputs={}, outputs={}, backup_comps={}
):

    return transform.model.ini_mapper(
        transform_type, inputs, outputs, backup_comps
    )
