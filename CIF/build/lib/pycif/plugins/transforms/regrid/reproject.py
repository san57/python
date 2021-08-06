from osgeo import gdal, ogr, osr

import numpy as np
from pyproj import Proj, transform

from pycif.utils.check import info


def reproject_emissions(
    emis_orig,
    zlonc_orig,
    zlatc_orig,
    zlonc_target,
    zlatc_target,
    resol=10,
    option="mean",
    wk_proj="wgs84",
    orig_regular=True,
    return_weight=False,
):

    # Calculate the number of bands
    if len(emis_orig.shape) == 2:
        info("Reprojecting 2D data")
        nbands = 1
        emis_orig = emis_orig[np.newaxis, ...]

    elif len(emis_orig.shape) == 3:
        info("Reprojecting 3D data")
        nbands = emis_orig.shape[0]

    info("Number of bands:", nbands)

    # Shifting longitudes beyond 180 degree east
    if wk_proj == "wgs84":
        zlonc_orig = (zlonc_orig + 180) % 360 - 180
        zlonc_target = (zlonc_target + 180) % 360 - 180

    # Projection definition
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromProj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    wgs84_pyproj = Proj(init="epsg:4326")

    if wk_proj == "wgs84":
        ref_proj = osr.SpatialReference()
        ref_proj.ImportFromProj4(
            "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        )
        ref_pyproj = Proj(init="epsg:4326")

    else:
        ref_proj = osr.SpatialReference()
        ref_proj.ImportFromProj4(wk_proj)
        ref_pyproj = Proj(wk_proj)

    # Convert lon/lat domain coordinates to reference coordinate system
    xc_orig, yc_orig = transform(
        wgs84_pyproj, ref_pyproj, zlonc_orig, zlatc_orig
    )
    xc_target, yc_target = transform(
        wgs84_pyproj, ref_pyproj, zlonc_target, zlatc_target
    )

    # Domains dimensions
    nzonal_orig, nmerid_orig = zlonc_orig[:-1, :-1].shape
    nzonal_target, nmerid_target = zlonc_target[:-1, :-1].shape

    # Rasterize to regular grid the origin emissions if necessary
    if not orig_regular:
        # Rasterize to regular grid the origin emissions if necessary
        (
            vector_grid_driver,
            vector_grid_ds,
            vector_grid_layer,
        ) = domain2polygons(zlonc_orig, zlatc_orig, wk_proj, data=emis_orig)

        # Rasterize irregular origin data
        x_min, x_max, y_min, y_max = vector_grid_layer.GetExtent()
        pixel_size = (
            max(
                np.median(np.ediff1d(xc_target)),
                np.median(np.ediff1d(xc_target.T)),
            )
            / resol
        )
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)

        geotransform = (x_min, pixel_size, 0, y_max, 0, -pixel_size)

        output_raster = gdal.GetDriverByName("MEM").Create(
            "", x_res, y_res, nbands, gdal.GDT_Float32
        )
        output_raster.SetGeoTransform(geotransform)
        output_raster.SetProjection(ref_proj.ExportToWkt())

        for iband in range(nbands):
            gdal.RasterizeLayer(
                output_raster,
                [iband + 1],
                vector_grid_layer,
                None,
                None,
                [0],
                options=["ATTRIBUTE=field_{}".format(iband)],
            )

    else:
        xmin, ymin, xmax, ymax = [
            xc_orig.min(),
            yc_orig.min(),
            xc_orig.max(),
            yc_orig.max(),
        ]
        pixel_xsize = max(
            np.median(np.ediff1d(xc_orig)), np.median(np.ediff1d(xc_orig.T))
        )
        pixel_ysize = max(
            np.median(np.ediff1d(yc_orig)), np.median(np.ediff1d(yc_orig.T))
        )
        x_npixl = int((xmax - xmin) / pixel_xsize)
        y_npixl = int((ymax - ymin) / pixel_ysize)

        # Transposing lon/lat if not in expected order
        lonlat_transpose = False
        if np.abs(x_npixl - nzonal_orig) > np.abs(x_npixl - nmerid_orig):
            lonlat_transpose = True
            emis_orig = np.transpose(emis_orig, axes=(0, 2, 1))
            xc_orig = xc_orig.T
            yc_orig = yc_orig.T
            nzonal_orig, nmerid_orig = zlonc_orig[:-1, :-1].shape

        # Shift original array spanning beyond 180 degrees east
        xpxl_ref = None
        if wk_proj == "wgs84":
            if (
                xc_orig[0, 0] + nzonal_orig * (xc_orig[1, 1] - xc_orig[0, 0])
                > 180
            ):
                xpxl_ref = int(
                    np.round(
                        (180 - xc_orig[0, 0]) / (xc_orig[1, 1] - xc_orig[0, 0])
                    )
                )

                xc_orig = np.concatenate(
                    (xc_orig[xpxl_ref:], xc_orig[:xpxl_ref]), axis=0
                )
                yc_orig = np.concatenate(
                    (yc_orig[xpxl_ref:], yc_orig[:xpxl_ref]), axis=0
                )
                emis_orig = np.concatenate(
                    (emis_orig[:, xpxl_ref:], emis_orig[:, :xpxl_ref]), axis=1
                )

        geotransform = (
            xc_orig[0, 0],
            xc_orig[1, 1] - xc_orig[0, 0],
            0,
            yc_orig[0, 0],
            0,
            yc_orig[1, 1] - yc_orig[0, 0],
        )

        # Create Raster with all emissions bands
        output_raster = gdal.GetDriverByName("MEM").Create(
            "", nzonal_orig, nmerid_orig, nbands, gdal.GDT_Float32
        )  #

        output_raster.SetGeoTransform(geotransform)

        # Exports the coordinate system to the file
        output_raster.SetProjection(ref_proj.ExportToWkt())

        # Loop over month to fill the raster
        for iband in range(nbands):
            # Writes my array to the raster
            output_raster.GetRasterBand(iband + 1).WriteArray(emis_orig[iband])

    # Create polygons from target domain
    vector_grid_driver, vector_grid_ds, vector_grid_layer = domain2polygons(
        zlonc_target, zlatc_target, wk_proj
    )

    # Compute projection
    emis_target = loop_zonal_stats(
        vector_grid_layer,
        output_raster,
        resol=resol,
        option=option,
        return_weight=return_weight,
    )

    # Return only weights if required
    if return_weight:
        if xpxl_ref is not None:
            for k, wgt in enumerate(emis_target):
                mask = wgt[1] >= xpxl_ref
                emis_target[k][1][mask] -= xpxl_ref
                emis_target[k][1][~mask] += xpxl_ref

        if lonlat_transpose:
            for k, wgt in enumerate(emis_target):
                emis_target[k] = (
                    emis_target[k][1],
                    emis_target[k][0],
                    emis_target[k][2],
                )

        return emis_target

    emis_target = np.reshape(emis_target, (nzonal_target, nmerid_target, -1))
    emis_target = np.transpose(emis_target, axes=(2, 1, 0))

    # Removes NaNs
    emis_target[np.isnan(emis_target)] = 0.0

    # if only one band, returns only the 2D dataset,
    # otherwise, returns everything
    if nbands == 1:
        return emis_target[0]
    else:
        return emis_target


def domain2polygons(zlonc_target, zlatc_target, wk_proj="wgs84", data=None):

    nzonal_target, nmerid_target = zlonc_target[:-1, :-1].shape
    if data is not None:
        nbands = data.shape[0]

    # Reference projections
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromProj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    wgs84_pyproj = Proj(init="epsg:4326")

    if wk_proj == "wgs84":
        ref_proj = osr.SpatialReference()
        ref_proj.ImportFromProj4(
            "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        )
        ref_pyproj = Proj(init="epsg:4326")

    else:
        ref_proj = osr.SpatialReference()
        ref_proj.ImportFromProj4(wk_proj)
        ref_pyproj = Proj(wk_proj)

    # Convert lon/lat domain coordinates to reference coordinate system
    xc, yc = transform(wgs84_pyproj, ref_pyproj, zlonc_target, zlatc_target)

    # define the layer
    vector_grid_driver = ogr.GetDriverByName("MEMORY")
    vector_grid_ds = vector_grid_driver.CreateDataSource("memData")
    vector_grid_layer = vector_grid_ds.CreateLayer(
        "grid", ref_proj, geom_type=ogr.wkbPolygon
    )

    if data is not None:
        for iband in range(nbands):
            idField = ogr.FieldDefn("field_{}".format(iband), ogr.OFTReal)
            vector_grid_layer.CreateField(idField)

    featureDefn = vector_grid_layer.GetLayerDefn()

    for i in range(nmerid_target):
        for j in range(nzonal_target):
            lon_tmp = [
                xc[j, i],
                xc[j, i + 1],
                xc[j + 1, i + 1],
                xc[j + 1, i],
                xc[j, i],
            ]

            if wk_proj == "wgs84":
                lon_tmp = np.unwrap(lon_tmp, discont=180)

            # Create ring
            ring = ogr.Geometry(ogr.wkbLinearRing)
            _ = ring.AddPoint(lon_tmp[0], yc[j, i])
            _ = ring.AddPoint(lon_tmp[1], yc[j, i + 1])
            _ = ring.AddPoint(lon_tmp[2], yc[j + 1, i + 1])
            _ = ring.AddPoint(lon_tmp[3], yc[j + 1, i])
            _ = ring.AddPoint(lon_tmp[4], yc[j, i])

            # Create polygon
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # Add feature to layer
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(poly)

            # Add data if any
            if data is not None:
                for iband in range(nbands):
                    feature.SetField(
                        "field_{}".format(iband), np.double(data[iband, j, i])
                    )

            vector_grid_layer.CreateFeature(feature)
            feature = None

    return vector_grid_driver, vector_grid_ds, vector_grid_layer


####################
# ZONAL STATISTICS #
####################
def zonal_stats(feat, raster, resol=10, option="mean", return_weight=False):
    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())

    # Create temporary layer from feature
    drv = ogr.GetDriverByName("MEMORY")
    ds = drv.CreateDataSource("memData")
    layer = ds.CreateLayer("grid", raster_srs, geom_type=ogr.wkbPolygon)
    layer.CreateFeature(feat)

    # Get raster georeference info
    geotransform = raster.GetGeoTransform()
    xOrigin = geotransform[0]
    yOrigin = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    xWest = min(xOrigin, xOrigin + pixelWidth * raster.RasterYSize)
    xEast = max(xOrigin, xOrigin + pixelWidth * raster.RasterYSize)
    ySouth = min(yOrigin, yOrigin + pixelHeight * raster.RasterXSize)
    yNorth = max(yOrigin, yOrigin + pixelHeight * raster.RasterXSize)

    # Get extent of feat
    geom = feat.GetGeometryRef()
    if geom.GetGeometryName() == "MULTIPOLYGON":
        count = 0
        pointsX = []
        pointsY = []
        for _ in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)
            count += 1
    elif geom.GetGeometryName() == "POLYGON":
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []
        pointsY = []
        for p in range(numpoints):
            lon, lat, z = ring.GetPoint(p)
            pointsX.append(lon)
            pointsY.append(lat)

    else:
        raise Exception(
            "ERROR: Geometry needs to be either Polygon or Multipolygon"
        )

    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    if xmax < xWest or xmin > xEast or ymax < ySouth or ymin > yNorth:
        info("Outside input domain")
        info(geotransform)
        info(pointsX)
        info(pointsY)
        return raster.RasterCount * [np.nan]

    # Specify offset and rows and columns to read
    xoff1 = int(np.floor((xmin - xOrigin) / pixelWidth))
    xoff2 = int(np.floor((xmax - xOrigin) / pixelWidth))
    xoff = max(0, min(xoff1, xoff2))

    yoff1 = int(np.floor((ymin - yOrigin) / pixelHeight))
    yoff2 = int(np.floor((ymax - yOrigin) / pixelHeight))
    yoff = max(0, min(yoff1, yoff2))

    xcount = min(max(xoff1, xoff2) + 1, raster.RasterYSize - 1) - xoff + 2
    ycount = min(max(yoff1, yoff2) + 1, raster.RasterXSize - 1) - yoff + 2

    # Create memory target raster
    target_ds = gdal.GetDriverByName("MEM").Create(
        "", ycount * resol, xcount * resol, 1, gdal.GDT_Byte
    )
    target_ds.SetGeoTransform(
        (
            xOrigin + xoff * pixelWidth,
            pixelWidth / resol,
            0,
            yOrigin + yoff * pixelHeight,
            0,
            pixelHeight / resol,
        )
    )

    # Rasterize zone polygon to raster
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])
    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray().astype(np.float)
    datamask2 = [
        [datamask[i::resol, j::resol] for i in range(resol)]
        for j in range(resol)
    ]
    datamask2 = np.array(datamask2).mean(axis=(0, 1))

    # If asks weights only
    if return_weight:
        meshi, meshj = np.meshgrid(
            np.arange(yoff, yoff + ycount), np.arange(xoff, xoff + xcount)
        )
        mask = datamask2 > 0

        return (meshi[mask], meshj[mask], datamask2[mask] / datamask2.sum())

    # Equivalent clipped raster
    nbands = raster.RasterCount

    raster_locmean = []
    for iband in range(1, nbands + 1):
        banddataraster = raster.GetRasterBand(iband)
        dataraster = banddataraster.ReadAsArray(
            yoff, xoff, ycount, xcount
        ).astype(np.float)

        # Averaging over the area corresponding to the pixel
        zoneraster = dataraster * datamask2 / datamask2.sum()
        if option == "mean":
            raster_locmean.append(np.sum(zoneraster))
        elif option == "median":
            raster_locmean.append(np.median(zoneraster))
        elif option == "min":
            raster_locmean.append(zoneraster[np.argmin(datamask2)])
        elif option == "max":
            raster_locmean.append(
                dataraster[
                    np.unravel_index(datamask2.argmax(), datamask2.shape)
                ]
            )
        else:
            raster_locmean.append(np.nan)

    # Calculate statistics of zonal raster
    return raster_locmean


def loop_zonal_stats(
    lyr, raster, resol=10, option="mean", return_weight=False
):
    featList = range(lyr.GetFeatureCount())
    statList = []

    for FID in featList:
        try:
            feat = lyr.GetFeature(FID)
            meanValue = zonal_stats(
                feat,
                raster,
                resol=resol,
                option=option,
                return_weight=return_weight,
            )
        except Exception as e:
            info(e)
            meanValue = raster.RasterCount * [np.nan]
            raise e

        if FID % (len(featList) / 10) == 0:
            info(FID, len(featList))

        statList.append(meanValue)

    if return_weight:
        return statList

    else:
        return np.array(statList)
