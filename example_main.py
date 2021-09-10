######################################################################
# Author: Christoph Schaller, BFH-HAFL, December 2020
#
# Script for demonstrating the use of pyFINT. The script expects
# input rasters of 1m resolution
######################################################################

import os
import sys
import importlib
import argparse

from datetime import datetime, date, time, timedelta
import time
from pyfintcontroller import *

# Default entry point
if __name__ == "__main__":

    start_time = time.time()

    #Expected resolution is 1m
    #One Filter example
    #One 1.5m resize example

    #Change to select, which test is run
    variant = 3

    #Path to output folder
    working_dir = os.getcwd()
    #Paths to input rasters 
    #Vegetation Height Model/Normalised Surface Model
    nsm_file = "VHM_1m.tif" 

    #
    # Standard Detection with 1m input VHM without resizing of filtering
    #
    fint_controller = pyFintController()
    fint_controller.set_working_dir(working_dir)
    #Whether to allow the use of altitude in DBH calculation (requires DTM)
    fint_controller.m_altitude_allowed = False
    #NSM/VHM used for detection
    fint_controller.set_normalized_model_file_name(nsm_file,None)
    #Tell the controller to use the NSM as source
    fint_controller.use_normalized_surface_model_as_input(True)
    #Set the function for calculating the DBH, whether to allow altitude in calculation
    fint_controller.set_dbh_function("2.52*H^0.84", False)
    #Whether to randomize the DBH value and the degree of deviation in percent
    fint_controller.set_diameter_randomization(False,20) 
    #Minimum height of a pixel to be considered for a local maxima
    fint_controller.set_minimum_height(1)
    #Minimum height for a detected maxima to be consideres as a tree
    fint_controller.set_minimum_detection_height(4)
    #Tell the controller to run the detection
    fint_controller.run_process()

    #
    # Detection with 1m input VHM resized to 1.5m
    #
    fint_controller = pyFintController()
    fint_controller.set_working_dir(working_dir)
    #Whether to allow the use of altitude in DBH calculation (requires DEM)
    fint_controller.m_altitude_allowed = False
    #NSM/VHM used for detection
    fint_controller.set_normalized_model_file_name(nsm_file,None)
    #Tell the controller to use the NSM as source
    fint_controller.use_normalized_surface_model_as_input(True)
    #Set the function for calculating the DBH, whether to allow altitude in calculation
    fint_controller.set_dbh_function("2.52*H^0.84", False)
    #Whether to randomize the DBH value and the degree of deviation in percent
    fint_controller.set_diameter_randomization(False,20) 
    #Minimum height of a pixel to be considered for a local maxima
    fint_controller.set_minimum_height(1)
    #Minimum height for a detected maxima to be consideres as a tree
    fint_controller.set_minimum_detection_height(4)
    #Tell the controller to resize the input tho the specified resolution with the given method
    #Supported methods basing on gdal: ["near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3"]
    fint_controller.set_resize_resolution(1.5,"bilinear")
    #Tell the controller to run the detection
    fint_controller.run_process()

    #
    # Detection with 1m input VHM and with Gauss filter sigma=2 and radius=3
    #
    fint_controller = pyFintController()
    fint_controller.set_working_dir(working_dir)
    #Whether to allow the use of altitude in DBH calculation (requires DEM)
    fint_controller.m_altitude_allowed = False
    #NSM/VHM used for detection
    fint_controller.set_normalized_model_file_name(nsm_file,None)
    #Tell the controller to use the NSM as source
    fint_controller.use_normalized_surface_model_as_input(True)
    #Set the function for calculating the DBH, whether to allow altitude in calculation
    fint_controller.set_dbh_function("2.52*H^0.84", False)
    #Whether to randomize the DBH value and the degree of deviation in percent
    fint_controller.set_diameter_randomization(False,20) 
    #Minimum height of a pixel to be considered for a local maxima
    fint_controller.set_minimum_height(1)
    #Minimum height for a detected maxima to be consideres as a tree
    fint_controller.set_minimum_detection_height(4)
    #Tell the controller to apply a Gauss filter of the given strength and radius; radius needs to be an odd number
    fint_controller.set_gauss_filter(size = 3, sigma = 2)
    #Tell the controller to run the detection
    fint_controller.run_process()

    #
    # Detection with 1m input VHM with resizing to 1.5 as well as with Gauss filter sigma=2 and radius=3
    #
    fint_controller = pyFintController()
    fint_controller.set_working_dir(working_dir)
    #Whether to allow the use of altitude in DBH calculation (requires DEM)
    fint_controller.m_altitude_allowed = False
    #NSM/VHM used for detection
    fint_controller.set_normalized_model_file_name(nsm_file,None)
    #Tell the controller to use the NSM as source
    fint_controller.use_normalized_surface_model_as_input(True)
    #Set the function for calculating the DBH, whether to allow altitude in calculation
    fint_controller.set_dbh_function("2.52*H^0.84", False)
    #Whether to randomize the DBH value and the degree of deviation in percent
    fint_controller.set_diameter_randomization(False,20) 
    #Minimum height of a pixel to be considered for a local maxima
    fint_controller.set_minimum_height(1)
    #Minimum height for a detected maxima to be consideres as a tree
    fint_controller.set_minimum_detection_height(4)
    #Tell the controller to resize the input tho the specified resolution with the given method
    #Supported methods basing on gdal: ["near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3"]
    fint_controller.set_resize_resolution(1.5,"bilinear")
    #Tell the controller to apply a Gauss filter of the given strength and radius; radius needs to be an odd number
    fint_controller.set_gauss_filter(size = 3, sigma = 2)
    #Tell the controller to run the detection
    fint_controller.run_process()



    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))
