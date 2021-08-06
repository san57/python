from pycif.utils.path import init_dir
import os
from shutil import copytree, ignore_patterns, rmtree, copy


def ini_mapper(model, transform_type, inputs={}, outputs={}, backup_comps={}):

    default_dict = {'input_dates': model.input_dates, 'force_read': True,
                    'force_dump': True}
    dict_surface = dict(default_dict, **{'domain': model.domain})

    # Executable
    mapper = {'inputs':
                  {('fluxes', s): dict_surface
                   for s in ['CH4']},
              'outputs': {('concs', s): {}
                          for s in ['CH4']}
              }

    return mapper

