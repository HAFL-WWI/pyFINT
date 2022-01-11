######################################################################
# Copyright (C) 2021 BFH
#
# Helper functions for processing with pyFINT. 
#
# Author: Christoph Schaller, BFH-HAFL, December 2021
######################################################################


import os
import sys
import math

import rasterio
from rasterio import windows
import rasterio.mask
import rasterio.merge 

from osgeo import ogr, gdal

import logging

class LogFile(object):
    """File-like object to log text using the `logging` module."""
    def __init__(self, name=None):
        self.logger = logging.getLogger(name)

    def write(self, msg, level=logging.INFO):
        self.logger.log(level, msg)

    def flush(self):
        for handler in self.logger.handlers:
            handler.flush()


# Function to ensure that a directory exists
def ensure_dir (path):
    if not os.path.isdir(path):
        os.mkdir(path)

# Function that returns a list of immediate subdirectories
def get_subdirs (path):
    dirs = []
    for f in os.listdir(path):
        full_path = os.path.join(path,f)
        if(os.path.isdir(full_path)):
          dirs.append([f,full_path])

    return dirs

# Function that deletes an existing Shapefile
#Based on https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#delete-a-file
def delete_shapefile(path):
    DriverName = "ESRI Shapefile"
    driver = ogr.GetDriverByName(DriverName)
    if os.path.exists(path):
        driver.DeleteDataSource(path)

# Function that deletes an existing Geotiff
def delete_raster(raster):
    data = gdal.Open(raster, gdal.GA_ReadOnly)
    driver = data.GetDriver()
    data = None
    if os.path.exists(raster):
        driver.Delete(raster)

#Function to copy a GeoTIFF raster
def copy_raster_tiff(in_raster, out_raster):
    driver = gdal.GetDriverByName('GTiff')
    in_ds = gdal.Open(in_raster)
    out_ds = driver.CreateCopy(out_raster, in_ds, 0)
    in_ds = None
    out_ds = None

# Function for cropping a raster based on a set of vector geometries
def crop_image(input_path, output_path, mask, dtype=None):
    with rasterio.open(input_path,mode='r') as input_image:
        no_data = input_image.nodata
        out_image, out_transform = rasterio.mask.mask(input_image, mask,
                                                        all_touched=True, nodata=no_data, crop=True)


        out_meta = input_image.meta.copy()
        out_meta.update({#"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform,
                         "compress":"lzw"})
        if dtype:
            out_meta.update({"dtype":dtype})
        if no_data:
            out_meta.update({"nodata":no_data})

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)

#Function that returns the 2D euclidean distance  between two trees
def dist_2d(tree1, tree2):
    dx = tree1["X_LV95"]-tree2["X_LV95"]
    dy = tree1["Y_LV95"]-tree2["Y_LV95"]
    d = math.sqrt(dx**2+dy**2)
    return d
