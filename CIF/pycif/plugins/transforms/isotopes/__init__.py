from .native2obsvect import native2obsvect
from .native2state import native2state
from .obsvect2native import obsvect2native
from .state2native import state2native


def ini_mapper(transform, mapper, transform_type):
    """ Re-arrange the mapper depending on the isotopic transforms to carry out


    :param self:
    :param transform:
    :param mapper:
    :param transform_type:
    :return:
    """

    # Default names for signatures corresponding to different types of inputs
    if transform_type == "fluxes":
        input_signature = "isosignatures"

    else:
        input_signature = "iso{}".format(transform_type)

    # Loop over isotopic transforms
    for isostep in transform.isospecies.attributes:
        isostp = getattr(transform.isospecies, isostep)
        names = isostp.inputs.names
        outputs = isostp.outputs.names
        inputs = [input_signature, transform_type]

        name_iso = names[0]
        name_flx = names[1]

        if name_flx in mapper["output_parameter"]:
            ind = mapper["output_parameter"].index(name_flx)
            del mapper["output_parameter"][ind]
            del mapper["output_component"][ind]

        for (comp, param) in zip(inputs, names):
            if (comp, param) not in list(
                zip(mapper["input_component"], mapper["input_parameter"])
            ):
                mapper["input_component"].append(comp)
                mapper["input_parameter"].append(param)

        mapper["output_parameter"].extend(outputs)
        mapper["output_component"].extend(2 * [transform_type])

        # Temporary fix for other parameters
        mapper["force_read"].extend(2 * [True])
        mapper["input_dates_parameters"].extend(2 * [0])
        mapper["input_files_parameters"].extend(2 * [0])

        isostp.mapper = {
            "input_parameter": [name_iso, name_flx],
            "input_component": [input_signature, transform_type],
            "output_component": 2 * [transform_type],
            "output_parameter": outputs,
        }

    return mapper
