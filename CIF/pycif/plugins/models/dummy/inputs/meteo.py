def make_meteo(self, datastore, ddi, ddf):

    ds = datastore.datastore[("meteo", "")]["spec"]
    self.meteo.data = ds.loc[(ds.index >= ddi) & (ds.index <= ddf)]
