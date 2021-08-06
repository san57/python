from types import MethodType

from .simul import simul

requirements = {'controlvect': {'any': True, 'empty': False},
                'obsvect': {'any': True, 'empty': False},
                'obsoperator': {'any': True, 'empty': True,
                                'name': 'standard', 'version': 'std'}
                }
