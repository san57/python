from pycif.utils.check import verbose


def check_inputs(inputs, mode):
    """Check the consistency of inputs given to the observation operator.
    
    """
    if mode not in ['tl', 'fwd', 'adj']:
        verbose("The following running mode is not accepted by the "
                "observation operator: {}".format(mode))
        raise Exception
    
    if mode == 'tl' and not (hasattr(inputs, 'x') and hasattr(inputs, 'dx')):
        verbose("The observation operator was operated in tangent-linear mode "
                "but not with both increments and control vector")
        raise Exception
    
    if mode == 'fwd' and not hasattr(inputs, 'x'):
        verbose("The observation operator was operated in forward mode "
                "with no control vector")
        verbose("All inputs will be dealt as fixed")
    
    return True
