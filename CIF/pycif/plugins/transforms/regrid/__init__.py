import os
from shutil import copy

from pycif.utils.path import init_dir
from .state2native import state2native

requirements = {"domain": {"any": True, "empty": False}}


def ini_data(plugin, **kwargs):

    dir_regrid = "{}/transforms/regrid/".format(plugin.workdir)
    plugin.dir_regrid = "{}/transforms/regrid/".format(plugin.workdir)

    # Initializes reference directories if needed
    init_dir(dir_regrid)


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

    # Copy weigth file if specified
    for id_in in mapper["inputs"]:
        basename = "wgt_file_{}_{}.pickle".format(*id_in)
        target_file = "{}/{}".format(transform.dir_regrid, basename)
        orig_file = "{}/{}".format(getattr(transform, "dir_wgt", ""), basename)
        if os.path.isfile(orig_file):
            copy(orig_file, target_file)
        mapper["inputs"][id_in]["weight_file"] = target_file

    return mapper
