######################################################################
# Copyright (C) 2021 BFH
#
# Script for processing Cantons in the FINT-CH project. 
#
# Author: Christoph Schaller, BFH-HAFL, December 2021
######################################################################

import os
import sys
import math
import glob

import numpy as np

import fiona

from shapely.geometry import Point, box
from osgeo import ogr

#Path to the folder containing pyFINT
PYFINT_HOME = os.environ.get("PYFINT_HOME")
sys.path.append(PYFINT_HOME)
from pyfintcontroller import *

#Path to the folder containing the FINT-CH artifacts
FINTCH_HOME = os.environ.get("FINTCH_HOME")
sys.path.append(os.path.join(FINTCH_HOME,"Common"))
from fintch_utilities import *

import numpy as np
import pandas as pd
import geopandas as gpd

from osgeo import ogr, osr, gdal
import psycopg2
from shapely.wkt import dumps, loads

import configparser

from datetime import datetime, date, time, timedelta
import time

from multiprocessing import Process, Pool, Queue, JoinableQueue, current_process, freeze_support

import logging
import traceback

import fintch_processing_core

def worker(q, work_function, cfg):
    db_connection = psycopg2.connect(host=cfg.get("db","host"), dbname=cfg.get("db","dbname"), user=cfg.get("db","user"), password=cfg.get("db","password"))

    configure_log(cfg)
    
    current_forest_mask_path = None
    current_trasse_mask_path = None
    forest_mask_df = None

    while True:
        #Consume work as long as there are items in the queue
        if q.qsize()>0:
            flaeche_record = q.get()
            if "waldmaske" in flaeche_record and flaeche_record["waldmaske"] != None and flaeche_record["waldmaske"] != "" :
                if flaeche_record["waldmaske"] != current_forest_mask_path:
                    current_forest_mask_path = flaeche_record["waldmaske"]
                    forest_mask_df = gpd.read_file(current_forest_mask_path)

                flaeche_record["waldmaske_df"] = forest_mask_df

            if "trasse_maske" in flaeche_record and flaeche_record["trasse_maske"] != None and flaeche_record["trasse_maske"] != "":
                if flaeche_record["trasse_maske"] != current_trasse_mask_path:
                    current_trasse_mask_path = flaeche_record["trasse_maske"]
                    trasse_mask_df = gpd.read_file(current_trasse_mask_path)

                flaeche_record["trasse_mask_df"] = trasse_mask_df


            work_function(flaeche_record,db_connection)
            q.task_done()

        else:
            #No more work available
            print("Exit:",current_process())
            db_connection.close()
            return

def process_setup(parameter_sets, flaechen_df, flaeche_id_column, flaeche_info_df, process_function, table_schema, table_base_name, cfg, add_forest_mask = False, add_trasse_mask = False, num_processes = 1):
    # Create queues
    perimeter_queue = JoinableQueue()

    #Execute detection by perimeter/flaeche
    for i, flaeche in flaeche_info_df.iterrows():
        
        vhm_path = flaeche["VHM"]
        vhm_input_file = vhm_path

        flaeche_id = flaeche["Flaeche_ID"]
        flaeche_geom = flaechen_df[flaechen_df[flaeche_id_column]==flaeche["Flaeche_ID"]].geometry.unary_union #A perimneter can have several parts
        flaeche_bounds = flaeche_geom.bounds
        minx = int((flaeche_bounds[0]//1000)*1000)
        miny = int((flaeche_bounds[1]//1000)*1000)
        maxx = int((flaeche_bounds[2]//1000)*1000+1000)
        maxy = int((flaeche_bounds[3]//1000)*1000+1000)

        for x in range(minx,maxx,1000):
            for y in range(miny,maxy,1000):
                tile_id = x//1000*10000+y//1000

                tile_geom = box(x,y,x+1000,y+1000)
                if tile_geom.intersects(flaeche_geom): #Actually belongs to this perimeter

                    perimeter_id = tile_id
                    perimeter_record = {
                        "parameter_sets": parameter_sets,
                        "perimeter_id": perimeter_id,
                        "geom_flaeche":flaeche_geom,
                        "vhm_input_file": vhm_path,
                        "geometry": tile_geom,
                        "flaeche_id": flaeche_id,
                        "quelle_id": 0,
                        "result_base_path": cfg.get("fintch_processing_paths","result_base_path"),

                        "perimeter_buffer" : cfg.getint("pyfint","perimeter_buffer"),
                        "r_max" : cfg.getfloat("pyfint","r_max"),
                        "epsg" : cfg.get("pyfint","epsg"),
                        "crs" : {'init': "epsg:"+cfg.get("pyfint","epsg")},

                        "table_schema": table_schema, 
                        "table_base_name": table_base_name,
                        
                        "mischungsgrad": flaeche["Mischungsgrad"],
                        "vhm_input_file_150": flaeche["VHM_150"]
                    }
                    if add_forest_mask:
                        perimeter_record["waldmaske"] = flaeche["Waldmaske"]
                    if add_trasse_mask:
                        perimeter_record["trasse_maske"] = flaeche["Trasse"]

                    perimeter_queue.put(perimeter_record)

    #Create and start worker processes 
    for i in range(num_processes):
        proc = Process(target=worker, args=(perimeter_queue,process_function,cfg,))
        print("Start: ",proc)
        proc.start()

    perimeter_queue.join()

    print("Processing finished")

def configure_log(cfg):
    log_path = cfg.get("fintch_processing_paths","log_path")
    logfile_info_path = os.path.join(log_path, current_process().name+"_info.log")
    logfile_error_path = os.path.join(log_path, current_process().name+"_error.log")
    log_format = "%(asctime)s; %(processName)s; %(levelname)s; %(name)s; %(message)s"

    # comment this to suppress console output
    stream_handler = logging.StreamHandler()

    file_handler_info = logging.FileHandler(logfile_info_path, mode='w')
    file_handler_info.setLevel(logging.INFO)

    file_handler_error = logging.FileHandler(logfile_error_path, mode='w')
    file_handler_error.setLevel(logging.ERROR)

    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            stream_handler, file_handler_info, file_handler_error
        ])
    sys.stdout = LogFile('stdout')
    sys.stderr = LogFile('stderr')



# Default entry point
if __name__ == "__main__":

    start_time = time.time()

    #Setup detection
    path_to_config_file = os.environ['FINTCH_CONFIG_HOME']
    ini_config_file = os.path.join(path_to_config_file, "FINTCH_config.ini")

    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(ini_config_file)

    result_base_path = cfg.get("fintch_processing_paths","result_base_path")
    log_path = cfg.get("fintch_processing_paths","log_path")
    flaechen_info_path = r"Q:\fint-ch\Geodaten\kantone_info.csv"
    flaechen_path = r"Q:\fint-ch\Geodaten\Kantone\outline.shp"  
    flaeche_id_column = "fid"

    flaechen_df = gpd.read_file(flaechen_path)
    flaeche_info_df = pd.read_csv(flaechen_info_path, delimiter=";")

    ensure_dir(result_base_path)
    ensure_dir(log_path)

    truncate = False

    configure_log(cfg)

    table_schema = "fintch"
    table_base_name = "fint"
    table_owner = "geoserver"
    add_trasse_mask = True
    db_connection = psycopg2.connect(host=cfg.get("db","host"), dbname=cfg.get("db","dbname"), user=cfg.get("db","user"), password=cfg.get("db","password"))
    srid = cfg.get("pyfint","epsg")
    fintch_processing_core.create_db_tables(table_schema,table_base_name,table_owner,srid,db_connection)
    
    parameter_sets = {
        31: {"vhm_source":"VHM_ALS", "dbh_function":"2.52*H^0.84", "randomized":False, "random_variance":0, "altitutde_allowed":False, "use_normalized_surfacemodel":True, "minimum_detection_tree_height":1, "minimum_tree_height":4, "gauss_sigma":"",  "gauss_size":"", "resize_method":"bilinear", "resize_resolution":1, "output_suffix":"", "preprocessing":"", "postprocessing":""},
        32: {"vhm_source":"VHM_ALS", "dbh_function":"2.52*H^0.84", "randomized":False, "random_variance":0, "altitutde_allowed":False, "use_normalized_surfacemodel":True, "minimum_detection_tree_height":1, "minimum_tree_height":4, "gauss_sigma":"",  "gauss_size":"", "resize_method":"bilinear", "resize_resolution":1.5, "output_suffix":"", "preprocessing":"", "postprocessing":""},
        34: {"vhm_source":"VHM_ALS", "dbh_function":"2.52*H^0.84", "randomized":False, "random_variance":0, "altitutde_allowed":False, "use_normalized_surfacemodel":True, "minimum_detection_tree_height":1, "minimum_tree_height":4, "gauss_sigma":2,  "gauss_size":3, "resize_method":"bilinear", "resize_resolution":1, "output_suffix":"", "preprocessing":"", "postprocessing":""},
    }

    process_setup(parameter_sets,flaechen_df,flaeche_id_column,flaeche_info_df,fintch_processing_core.process_detection,table_schema,table_base_name,cfg, add_forest_mask=add_forest_mask, num_processes = 40)

    cursor = db_connection.cursor()
    cursor.execute("TRUNCATE TABLE "+table_schema+"."+table_base_name+"_processed_tree")
    db_connection.commit()
    cursor.close()

    process_setup(parameter_sets,flaechen_df,flaeche_id_column,flaeche_info_df,fintch_processing_core.process_fst,table_schema,table_base_name,cfg,add_forest_mask=add_forest_mask,num_processes = 15)
    process_setup(parameter_sets,flaechen_df,flaeche_id_column,flaeche_info_df,fintch_processing_core.process_trees,table_schema,table_base_name,cfg,add_forest_mask=add_forest_mask,add_trasse_mask=add_trasse_mask,num_processes = 40)

    db_connection.close()

    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))


