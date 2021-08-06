import numpy as np


def create_domain(domain,
                  **kwargs):
    """Creates a grid if needed

    Args:
        domain (dictionary): dictionary defining the domain.

    Returns:
         

    Notes:
        We assume that center coordinates are used in the domain definition 

    Todo:
        

    """
    
    nlon = domain.nlon
    nlat = domain.nlat
    
    dlon = (domain.xmax - domain.xmin) / nlon
    dlat = (domain.ymax - domain.ymin) / nlat
    
    # Corner coordinates
    lonc = np.arange(domain.xmin, domain.xmin + nlon * dlon, dlon)
    latc = np.arange(domain.ymin, domain.ymin + nlat * dlat, dlat)
    
    # Center coordinates
    lon = lonc + dlon / 2.
    lat = latc + dlat / 2.
    
    # Meshgrids
    zlon, zlat = np.meshgrid(lon, lat)
    
    lonc = np.append(lonc, lonc[-1] + dlon)
    latc = np.append(latc, latc[-1] + dlat)
    
    zlonc, zlatc = np.meshgrid(lonc, latc)
    
    domain.zlon_in = zlon
    domain.zlat_in = zlat
    domain.zlonc_in = zlonc
    domain.zlatc_in = zlatc
    domain.lon_in = lon
    domain.lat_in = lat
    
    # Adding global grid definition
    xmin_glob = domain.xmin_glob
    ymin_glob = domain.ymin_glob
    nlon_glob = domain.nlon_glob
    nlat_glob = domain.nlat_glob
    dlon_glob = domain.dx_glob
    dlat_glob = domain.dy_glob
    
    # Center coordinates
    lon_glob = np.arange(xmin_glob, xmin_glob + nlon_glob * dlon_glob,
                         dlon_glob) + dlon_glob / 2.
    lat_glob = np.arange(ymin_glob, ymin_glob + nlat_glob * dlat_glob,
                         dlat_glob) + dlat_glob / 2.
    zlon_glob, zlat_glob = np.meshgrid(lon_glob, lat_glob)
    
    # Corner coordinates
    lonc_glob = np.arange(xmin_glob, xmin_glob + nlon_glob * dlon_glob,
                          dlon_glob)
    latc_glob = np.arange(ymin_glob, ymin_glob + nlat_glob * dlat_glob,
                          dlat_glob)
    lonc_glob = np.append(lonc_glob, lonc_glob[-1] + dlon_glob)
    latc_glob = np.append(latc_glob, latc_glob[-1] + dlat_glob)
    
    zlonc_glob, zlatc_glob = np.meshgrid(lonc_glob, latc_glob)
    
    # Indices for inversion domain into global domain (relative to lower left
    #  corner)
    ix1 = np.argmin(np.abs(lonc_glob - domain.lon_in[0] - 1))
    iy1 = np.argmin(np.abs(latc_glob - domain.lat_in[0] - 1))
    ix2 = np.argmin(np.abs(lonc_glob - domain.lon_in[-1] - 1))
    iy2 = np.argmin(np.abs(latc_glob - domain.lat_in[-1] - 1))
    
    domain.ix1 = ix1
    domain.ix2 = ix2
    domain.iy1 = iy1
    domain.iy2 = iy2
    
    # Outside domain
    ilon, ilat = np.meshgrid(np.arange(nlon_glob), np.arange(nlat_glob))
    in_indexes = np.array([ilat[iy1:iy2, ix1:ix2], ilon[iy1:iy2, ix1:ix2]])
    raveled_indexes = \
        np.ravel_multi_index(in_indexes, zlat_glob.shape).flatten()
    domain.raveled_indexes_glob = raveled_indexes
    
    lat_out = np.delete(zlat_glob, raveled_indexes)
    lon_out = np.delete(zlon_glob, raveled_indexes)
    
    # Putting together longitudes and latitudes
    domain.zlon = np.append(domain.zlon_in.flatten(), lon_out)[:, np.newaxis]
    domain.zlat = np.append(domain.zlat_in.flatten(), lat_out)[:, np.newaxis]
    
    # Build all corners
    size_in = domain.zlon_in.size
    domain.zlatc = np.empty((4, domain.zlat.size))
    domain.zlonc = np.empty((4, domain.zlat.size))

    domain.zlatc[0, :size_in] = domain.zlat_in.flatten() - dlat / 2
    domain.zlatc[1, :size_in] = domain.zlat_in.flatten() - dlat / 2
    domain.zlatc[2, :size_in] = domain.zlat_in.flatten() + dlat / 2
    domain.zlatc[3, :size_in] = domain.zlat_in.flatten() + dlat / 2

    domain.zlonc[0, :size_in] = domain.zlon_in.flatten() - dlon / 2
    domain.zlonc[1, :size_in] = domain.zlon_in.flatten() + dlon / 2
    domain.zlonc[2, :size_in] = domain.zlon_in.flatten() + dlon / 2
    domain.zlonc[3, :size_in] = domain.zlon_in.flatten() - dlon / 2
    
    domain.zlatc[0, size_in:] = lat_out - dlat_glob / 2
    domain.zlatc[1, size_in:] = lat_out - dlat_glob / 2
    domain.zlatc[2, size_in:] = lat_out + dlat_glob / 2
    domain.zlatc[3, size_in:] = lat_out + dlat_glob / 2

    domain.zlonc[0, size_in:] = lon_out - dlon_glob / 2
    domain.zlonc[1, size_in:] = lon_out + dlon_glob / 2
    domain.zlonc[2, size_in:] = lon_out + dlon_glob / 2
    domain.zlonc[3, size_in:] = lon_out - dlon_glob / 2
    