from .state2native import state2native


def ini_mapper(
    transform, transform_type, inputs={}, outputs={}, backup_comps={}
):

    default_dict_in = {"force_read": True, "force_dump": False}
    default_dict_out = {"force_read": True}
    iter_ids = (
        [(transform.component[0], transform.parameter[0])]
        if len(transform.component) == 1
        else zip(transform.component, transform.parameter)
    )

    iter_dict_out = {
        trid: dict(inputs[trid], **default_dict_out) for trid in iter_ids
    }

    iter_dict_in = {
        trid: dict(inputs[trid], **default_dict_in) for trid in iter_ids
    }

    mapper = {
        "inputs": {trid: iter_dict_in[trid] for trid in iter_ids},
        "outputs": {trid: iter_dict_out[trid] for trid in iter_ids},
    }

    return mapper
