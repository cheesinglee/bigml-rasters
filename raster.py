# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 20:40:02 2014

@author: cheesinglee
"""

import os
import struct
import configparser
from numpy import array, linspace, meshgrid,reshape
import gdal
import osr

# dicts for constructing struct format string from R raster datatype
DATATYPES = {'LOG1S': '?', 'INT1S': 'b', 'INT2S': 'h',
             'INT4S': 'i', 'INT8S': 'q', 'INT1U': 'B',
             'INT2U': 'H', 'FLT4S': 'f', 'FLT8S': 'd'}
BYTEORDER = {'little': '<', 'big': '>'}

# lookup for GDAL datatypes from R raster datatype
GDALTYPES = {'LOG1S': gdal.GDT_Byte, 'INT1S': gdal.GDT_Byte,
             'INT2S': gdal.GDT_Int16, 'INT4S': gdal.GDT_Int32,
             'INT8S': gdal.GDT_Int32, 'INT1U': gdal.GDT_Byte,
             'INT2U': gdal.GDT_UInt16, 'FLT4S': gdal.GDT_Float32,
             'FLT8S': gdal.GDT_Float64}

# enable GDAL exceptions
gdal.UseExceptions()


def all_equal(sequence):
    return sequence.count(sequence[0]) == len(sequence)    


def load_r_raster(filename):
    config = configparser.ConfigParser()
    tokens = filename.split(os.path.extsep)

    tokens[-1] = 'grd'
    ini_path = os.path.extsep.join(tokens)

    tokens[-1] = 'gri'
    bin_path = os.path.extsep.join(tokens)

    res = config.read(ini_path)
    if len(res) == 0:
        raise(Exception, 'could not read R raster file %s' % ini_path)

    with open(bin_path, 'rb') as fid:
        binstring = fid.read()
    datatype = config.get('data', 'datatype')
    byteorder = config.get('data', 'byteorder')
    dt = DATATYPES[datatype]
    bo = BYTEORDER[byteorder]
    fmt = '%s%d%s' % (bo, len(binstring)/struct.calcsize(dt), dt)
    values = struct.unpack(fmt, binstring)

    ncols = config.getint('georeference', 'ncols')
    nrows = config.getint('georeference', 'nrows')
    xmax = config.getfloat('georeference', 'xmax')
    xmin = config.getfloat('georeference', 'xmin')
    ymax = config.getfloat('georeference', 'ymax')
    ymin = config.getfloat('georeference', 'ymin')

    xres = (xmax-xmin)/ncols
    yres = (ymax-ymin)/nrows

    minval = config.getfloat('data', 'minvalue')
    maxval = config.getfloat('data', 'maxvalue')
    nodatavalue = config.getfloat('data', 'nodatavalue')

    # make sure nodata values are correct
    values = [x if (x <= minval and x >= maxval) else nodatavalue
              for x in values]

    # pack into numpy array
    arr = array(values)
    arr = reshape(arr, (nrows, ncols))
    
    # create GDAL in-memory raster dataset
    driver = gdal.GetDriverByName('MEM')
    tokens[-1] = 'tif'
    raster_path = os.path.extsep.join(tokens)
    gdal_ds = driver.Create(raster_path, ncols, nrows, 1, GDALTYPES[datatype])
    gdal_ds.SetGeoTransform((xmin, xres, 0, ymax, 0, yres))
    band = gdal_ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.SetNoDataValue(nodatavalue)

    # projection
    proj4 = config.get('georeference', 'projection')
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj4)
    gdal_ds.SetProjection(srs.ExportToWkt())

    gdal_ds.FlushCache()    
        
    return gdal_ds


def load(filename, rastertype=None):
    if not rastertype:
        tokens = filename.split(os.path.extsep)
        ext = tokens[-1]
        if ext in ['gri', 'grd']:
            rastertype = 'r-raster'
        elif ext == 'bil':
            rastertype = 'bil'

    if rastertype == 'r-raster':
        return load_r_raster(filename)
    else:
        try:
            return gdal.Open(filename, gdal.GA_ReadOnly)
        except:
            print('could not open raster file %s' % filename)


def make_table(rasters):
    # get raster dimensions
    geotransform = rasters[0].GetGeoTransform()
    cols = rasters[0].RasterXSize
    rows = rasters[0].RasterYSize
    xres = geotransform[1] 
    yres = geotransform[5]
    topleft_x = geotransform[0]
    topleft_y = geotransform[3]
    
    instances = [dict() for i in range(cols*rows)]
    for r in rasters:
        # use the filename in description as feature name
        fn = r.GetDescription()
        name = os.path.basename(fn).split(os.path.extsep)[0]
        n_bands = r.RasterCount
        bands = [r.GetRasterBand(i+1) for i in range(n_bands)]
        for b in bands:
            name = b.GetDescription()
            if name == '':
                if n_bands > 1:
                    band_name = '%s_%d' % (name, i+1)
                else:
                    band_name = name
            values = b.ReadAsArray().flatten()
            for i, v in enumerate(values):
                instances[i][band_name] = v
    
    # add xy coordinates
    x_start = topleft_x + xres/2
    x_end = x_start + xres*cols
    y_start = topleft_y - yres/2
    y_end = y_start - yres*rows
    x_coords = linspace(x_start, x_end, cols)
    y_coords = linspace(y_start, y_end, rows)
    x_mesh, y_mesh = meshgrid(x_coords, y_coords)
    for i, inst in enumerate(instances):
        inst['x'] = x_mesh[i]
        inst['y'] = y_mesh[i]
    return instances

def rasters_to_table(filenames):
    rasters = [load(f) for f in filenames]
    geotransforms,projections = zip(*[(r.GetGeoTransform(),r.GetProjection()) 
                                    for r in rasters])
                                        
    if not all_equal(geotransforms) or not all_equal(projections):
        rasters = resample_and_reproject(rasters)
    return make_table(rasters)
        
if __name__ == '__main__':
    import glob
    filelist = glob.glob('data/bioclim/current/*.grd')
    rasters_to_table(filelist)
