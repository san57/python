from pycif.utils.check import info, debug
import copy


def native2inputs(self, datastore, input_type, datei, datef,
                  runsubdir, mode='fwd'):
    """Converts data at the model data resolution to model compatible input
    files.
    
    Args:
        self: the model Plugin
        input_type (str): one of 'fluxes', 'obs'
        datastore: data to convert
            if input_type == 'fluxes', a dictionary with flux maps
            if input_type == 'obs', a pandas dataframe with the observations
        datei, datef: date interval of the sub-simulation
        mode (str): running mode: one of 'fwd', 'adj' and 'tl'
        runsubdir (str): sub-directory for the current simulation
        workdir (str): the directory of the whole pycif simulation
    
    Notes:
        Copied from LMDZ. We do not attempt to run the model at this point.

        
    """
    
    ddi = min(datei, datef)
    ddf = max(datei, datef)
    
    if input_type not in ['fluxes', 'obs']:
        raise Exception("FLEXPART received unexpected input type: {}"
                        .format(input_type))
    
    elif input_type == 'obs':
        debug('Obs already stored in dataobs')
        return
    
    # Storing fluxes for later, copying to avoid erasing info elsewhere
    self.dataflx = {trid: {k: datastore.datastore[trid][k]
                           for k in datastore.datastore[trid]}
                    for trid in datastore.datastore}
    