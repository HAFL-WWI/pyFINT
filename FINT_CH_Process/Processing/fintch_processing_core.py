######################################################################
# Copyright (C) 2021 BFH
#
# Script with base routines for processing steps in the FINT-CH project. 
#
# Author: Christoph Schaller, BFH-HAFL, December 2021
######################################################################

import os
import sys
import math


from rasterstats import zonal_stats

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
from geopandas.tools import sjoin


from osgeo import ogr, osr, gdal
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *

import psycopg2
from shapely.wkt import dumps, loads
from shapely.geometry import Point, box


from datetime import datetime, date, time, timedelta
import time

from multiprocessing import Process, Pool, Queue, JoinableQueue, current_process, freeze_support

import logging
import traceback



def create_db_tables(table_schema, table_base_name, table_owner, srid, db_connection):
    """Method for creating the PostGIS database tables needed in the FINT-CH process.
    Existing tables are dropped beforehand.

    Args:
        table_schema (string): Name of the schema to create the tables in
        table_base_name (string): Base name for the created tables
        table_owner (string): Owner of the created tables
        db_connection (connection): psycopg2 connection to use for creating the tables
    """

    create_table_template = """
----
-- Table: raw detected trees
----
DROP TABLE IF EXISTS {0}.{1}_tree_detected;

CREATE TABLE {0}.{1}_tree_detected
(
    gid serial NOT NULL,
    x double precision,
    y double precision,
    hoehe real,
    dominanz real,
    bhd real,
    geom geometry(Point,{3}),
    parameterset_id smallint,
    perimeter_id integer,
    flaeche_id integer,
    hoehe_modified real,
    CONSTRAINT {1}_tree_detected_pkey PRIMARY KEY (gid)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE {0}.{1}_tree_detected
    OWNER to {2};

-- Index: geom
CREATE INDEX sidx_{1}_tree_detected_geom_idx
    ON {0}.{1}_tree_detected USING gist
    (geom)
    TABLESPACE pg_default;

-- Index parameterset_id
CREATE INDEX idx_{1}_tree_detected_parameterset_id
    ON {0}.{1}_tree_detected USING btree
    (parameterset_id)
    TABLESPACE pg_default;

-- Index parameterset_id, perimeter_id
CREATE INDEX idx_{1}_tree_detected_parameterset_id_perimeter_id
    ON {0}.{1}_tree_detected USING btree
    (parameterset_id, perimeter_id)
    TABLESPACE pg_default;

----
-- Table: detection perimeter
----
DROP TABLE IF EXISTS {0}.{1}_perimeter;

CREATE TABLE {0}.{1}_perimeter
(
    gid serial NOT NULL,
    geom geometry(Polygon,{3}),
    perimeter_id integer,
    flaeche_id integer,
    CONSTRAINT {1}_perimeter_pkey PRIMARY KEY (gid)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE {0}.{1}_perimeter
    OWNER to {2};

-- Index: geom
CREATE INDEX sidx_{1}_perimeter_geom
    ON {0}.{1}_perimeter USING gist
    (geom)
    TABLESPACE pg_default;
	
----	
-- Table: forest structure type raster
----
DROP TABLE IF EXISTS {0}.{1}_fst_raster;

CREATE TABLE {0}.{1}_fst_raster
(
    gid serial NOT NULL,
    geom geometry(Polygon,{3}),
    flaeche_id integer, 
    perimeter_id integer, 
    tile_id bigint,
    hdom smallint,
    dg smallint,
    nh smallint,
    fst smallint,
    CONSTRAINT {1}_fst_raster_pkey PRIMARY KEY (gid)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE {0}.{1}_fst_raster
    OWNER to {2};

-- Index: geom
CREATE INDEX sidx_{1}_fst_raster_geom_idx
    ON {0}.{1}_fst_raster USING gist
    (geom)
    TABLESPACE pg_default;

-- Index flaeche_id, perimeter_id
CREATE INDEX idx_{1}_fst_raster_flaeche_id_perimeter_id
    ON {0}.{1}_fst_raster USING btree
    (flaeche_id, perimeter_id)
    TABLESPACE pg_default;

----
-- Table: trees filtered by forest structure type
----
DROP TABLE IF EXISTS {0}.{1}_processed_tree;

CREATE TABLE {0}.{1}_processed_tree
(
    gid serial NOT NULL,
    x double precision,
    y double precision,
    hoehe real,
    dominanz real,
    bhd real,
    geom geometry(Point,{3}),
    parameterset_id smallint,
    fst_raster_id integer,
    flaeche_id integer,
    hoehe_modified real,
    fst smallint,
    CONSTRAINT {1}_processed_tree_pkey PRIMARY KEY (gid)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE {0}.{1}_processed_tree
    OWNER to {2};

-- Index: geom
CREATE INDEX sidx_{1}_processed_tree_geom_idx
    ON {0}.{1}_processed_tree USING gist
    (geom)
    TABLESPACE pg_default;
    """

    cursor = db_connection.cursor()
    sql = create_table_template.format(table_schema, table_base_name, table_owner, srid)
    cursor.execute(sql)
    db_connection.commit()
    cursor.close()




def generate_grid(min_x, min_y, max_x, max_y, out_shape_path, crs=2056, step_x=25,step_y=25):
    # create output file
    srs = osr.SpatialReference()
    srs.ImportFromEPSG( crs )

    logger = logging.getLogger()

    out_driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(out_shape_path):
        delete_shapefile(out_shape_path)
    out_ds = out_driver.CreateDataSource(out_shape_path)
    out_layer = None
    try:
        out_layer = out_ds.CreateLayer(out_shape_path,srs=srs,geom_type=ogr.wkbPolygon )
    except Error as ex: 
        logger.error("Error generating grid ", out_shape_path)
        logger.error(traceback.format_exception(*sys.exc_info()))
        raise ex


    out_layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64))
    feature_defn = out_layer.GetLayerDefn()

    cur_x = min_x
    cur_y = min_y
    col = 0
    row = -1

    e = len(str(int(min_x)))
    f = 10**e

    # create grid cells
    while cur_y < max_y:
        row += 1

        cur_x = min_x
        col = -1

        while cur_x < max_x:
            col += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(cur_x, cur_y)
            ring.AddPoint(cur_x+step_x, cur_y)
            ring.AddPoint(cur_x+step_x, cur_y+step_y)
            ring.AddPoint(cur_x, cur_y+step_y)
            ring.AddPoint(cur_x, cur_y)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            out_feature = ogr.Feature(feature_defn)
            out_feature.SetGeometry(poly)
            out_feature.SetField('id', cur_x*f+cur_y )
            out_layer.CreateFeature(out_feature)
            out_feature.Destroy

            cur_x += step_x

        cur_y += step_y

    # Close DataSources
    out_ds.Destroy()

def determine_fst(grid_path,vhm150_path,mixing_degree_path,envelope):

    output_folder = os.path.dirname(grid_path)

    logger = logging.getLogger()

    #Clip VHM
    vhm_output_file = os.path.join(output_folder,"vhm150_clip.tif")
    try:
        crop_image(vhm150_path, vhm_output_file, [envelope]) # Nodata Value may depend on source!
    except ValueError as ve:
        logger.error(grid_path)
        logger.error(traceback.format_exception(*sys.exc_info()))
        return -1
    
    #Clip NH
    mg_output_file = os.path.join(output_folder,"mg_clip.tif")
    try:
        crop_image(mixing_degree_path, mg_output_file, [envelope]) # Nodata Value may depend on source!
    except ValueError as ve:
        logger.error(grid_path)
        logger.error(traceback.format_exception(*sys.exc_info()))
        return -1

    ##
    ## Calculate hdom
    ##
    stats = zonal_stats(grid_path, vhm_output_file, stats=['percentile_80'], all_touched=True)

    # open grid polygon shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    grid_ds = driver.Open(grid_path, 1)
    layer = grid_ds.GetLayer()

    # add grid attribute fields
    layer.CreateField(ogr.FieldDefn('hdom', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('nh', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('dg_min', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('dg', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('FST', ogr.OFTInteger))

    # iterate over all features and add stand attribute values
    counter = 0
    for feature in layer:
        hp80 = 0
        if stats[counter].get('percentile_80') is not None:
            hp80 = stats[counter].get('percentile_80')

        # set and store hdom
        feature.SetField('hdom', hp80)
        layer.SetFeature(feature)
        counter += 1
    grid_ds, layer = None, None
    
    ##
    ## Calculate nh
    ##
    stats = zonal_stats(grid_path, mg_output_file, stats=['mean'], all_touched=True)

    grid_ds = driver.Open(grid_path, 1)
    layer = grid_ds.GetLayer()

    # iterate over all features and add stand attribute values
    counter = 0
    for feature in layer:
        nh = 0
        if stats[counter].get('mean') is not None:
            nh = stats[counter].get('mean')
            nh = round(nh/100)

        # set and store nh
        feature.SetField('nh', nh)
        layer.SetFeature(feature)
        counter += 1
    grid_ds, layer = None, None

    ##
    ## Calculate dg
    ##
    grid_ds = driver.Open(grid_path, 1)
    layer = grid_ds.GetLayer()

    # tmp files
    dg_classified_path = os.path.join(output_folder,"dg_layer.tif")
    tmp_lim_dg_path = os.path.join(output_folder,"dg_lim_dg.tif")

    # Layer threshold values (based on NFI definition, www.lfi.ch)
    min_height_hdom_factor_ms = 1.0 / 3.0
    min_height_hdom_factor_os = 2.0 / 3.0


    for feature in layer:
        # calculate and store dg_min
        hdom =  feature.GetFieldAsInteger('hdom')
        if hdom < 14: #Fix small stands issue
            dg_min = hdom*min_height_hdom_factor_ms
        else:
            dg_min = hdom*min_height_hdom_factor_os
        feature.SetField('dg_min', dg_min )
        layer.SetFeature(feature)
        counter += 1

    # Rasterize dg_min
    vhm_output_file

    vhm_ds = gdal.Open(vhm_output_file,GA_ReadOnly)

    driver_gtiff = gdal.GetDriverByName('GTiff')
    dg_min_ds = driver_gtiff.Create(tmp_lim_dg_path,vhm_ds.RasterXSize, vhm_ds.RasterYSize,1,gdal.GDT_Float32)
    dg_min_ds.SetGeoTransform(vhm_ds.GetGeoTransform())
    dg_min_ds.SetProjection(vhm_ds.GetProjection())

    dst_options = ['ATTRIBUTE=dg_min']
    gdal.RasterizeLayer(dg_min_ds, [1], layer, None, options=dst_options)

    # Produce "1" / "0" raster for each layer
    vhm_b1 = vhm_ds.GetRasterBand(1)
    dg_min_b1 = dg_min_ds.GetRasterBand(1)
    
    data_vhm = np.array(vhm_b1.ReadAsArray())
    data_dg_min = np.array(dg_min_b1.ReadAsArray())

    data_out = data_vhm>data_dg_min
    zoLembda = lambda x: 1 if x else 0
    vfunc = np.vectorize(zoLembda)
    data_out = vfunc(data_out)

    # Write the out file
    dst_options = ['COMPRESS=LZW']
    dg_ds = driver_gtiff.Create(dg_classified_path, vhm_ds.RasterXSize, vhm_ds.RasterYSize, 1, gdal.GDT_Byte, dst_options)
    CopyDatasetInfo(vhm_ds, dg_ds)
    band_out = dg_ds.GetRasterBand(1)
    BandWriteArray(band_out, data_out)

    vhm_ds, dg_min_ds, dg_ds = None, None, None

    # Zonal stats
    stats = zonal_stats(grid_path, dg_classified_path, stats=['mean'], all_touched=True)

    # iterate over all features and add stand attribute values
    counter = 0
    for feature in layer:
        dg = 0
        if stats[counter].get('mean') is not None:
            dg = stats[counter].get('mean')
            dg = round(dg*100)

        hdom = feature.GetFieldAsInteger('hdom')
        nh = feature.GetFieldAsInteger('nh')

        if nh <= 30:
            digit1 = 1
        elif 30 < nh <= 70:
            digit1 = 2
        else:
            digit1 = 3

        if dg <= 80:
            digit2 = 1
        else:
            digit2 = 2

        if hdom <= 22:
            digit3 = 1
        else:
            digit3 = 2

        fst = int(str(digit1) + str(digit2) + str(digit3))


        # set and store dg and FST
        feature.SetField('dg', dg)
        feature.SetField('FST', fst)
        layer.SetFeature(feature)
        counter += 1
    
    grid_ds, layer = None, None

    # Cleanup
    delete_raster(dg_classified_path)
    delete_raster(tmp_lim_dg_path)
    delete_raster(vhm_output_file)
    delete_raster(mg_output_file)
    
    return 0

def process_detection(record, db_connection):
    cursor = db_connection.cursor()

    fint_controller = pyFintController()
    logger = logging.getLogger()

    table_schema = record["table_schema"]
    table_base_name = record["table_base_name"]

    perimeter_insert_template = "INSERT INTO "+table_schema+"."+table_base_name+"_perimeter(geom, perimeter_id, flaeche_id) VALUES (ST_SetSRID(ST_GeomFromText('{0}'),{1}), {2}, {3});"
    tree_insert_template = "INSERT INTO "+table_schema+"."+table_base_name+"_tree_detected(x, y, hoehe, bhd, dominanz, geom, parameterset_id, perimeter_id, flaeche_id, hoehe_modified) VALUES ({0}, {1}, {2}, {3}, {4}, ST_SetSRID(ST_GeomFromText('{5}'),{6}), {7}, {8},{9},{10});" 

    result_base_path = record["result_base_path"]
    perimeter_buffer = record["perimeter_buffer"]
    r_max = record["r_max"]
    epsg = record["epsg"]
    crs = record["crs"]
    
    perimeter_id = record["perimeter_id"]
    flaeche_id = record["flaeche_id"]
    folder_name = "{0}_{1}".format(flaeche_id,perimeter_id)
    output_folder = os.path.join(result_base_path,folder_name)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    geom = record["geometry"]
    (minx, miny, maxx, maxy) = geom.bounds

    #Envelope by group
    bx = box(minx, miny, maxx, maxy)
    envelope = bx.buffer(perimeter_buffer, resolution=1).envelope

    sql = perimeter_insert_template.format(geom.wkt,epsg,perimeter_id,flaeche_id)
    cursor.execute(sql)
    db_connection.commit()

    parameter_sets = record["parameter_sets"]
    vhm_input_file = record["vhm_input_file"]

    fint_tree_dataframes = []

    for paramterset_id in parameter_sets:
        parameter_set = parameter_sets[paramterset_id]
        parameter_set["id"] = paramterset_id

        detection_result = detect_trees(parameter_set, output_folder, vhm_input_file, envelope, crs, fint_controller)
        if type(detection_result) == type(None):
            continue
        else:
            detection_result["parameterset_id"] = paramterset_id
            fint_tree_dataframes.append(detection_result)

    if len(fint_tree_dataframes)==0:
        cursor.close()
        return
    fint_trees_df = gpd.GeoDataFrame(pd.concat(fint_tree_dataframes, ignore_index=True), crs=fint_tree_dataframes[0].crs)

    geom_bounds = geom.bounds
    minx = geom_bounds[0]
    miny = geom_bounds[1]
    maxx = geom_bounds[2]
    maxy = geom_bounds[3]

    fint_trees_df = fint_trees_df[fint_trees_df["hoehe"].apply(lambda x: x.strip())!="nan"]
    fint_trees_df["hoehe"] = fint_trees_df["hoehe"].astype(np.double)
    fint_trees_df = fint_trees_df[(fint_trees_df["x"]>=minx) & (fint_trees_df["x"]<maxx) & (fint_trees_df["y"]>=miny) & (fint_trees_df["y"]<maxy) & (fint_trees_df.within(record["geom_flaeche"]))].copy()
    fint_trees_df["perimeter_id"] = perimeter_id
    fint_trees_df = merge_processing(fint_trees_df)

    if type(fint_trees_df)==type(None) or len(fint_trees_df) == 0:
        cursor.close()
        return
    fint_trees_df["wkt"] = fint_trees_df["geometry"].apply(dumps) #reflects onto the underlying df

    for tindex, ttree in fint_trees_df.iterrows():
        sql = tree_insert_template.format(ttree["x"],ttree["y"],ttree["hoehe"],ttree["bhd"],ttree["dominanz"],ttree["wkt"],epsg,ttree["parameterset_id"],perimeter_id,flaeche_id,ttree["hoehe_modified"] if ttree["hoehe_modified"].strip() != "nan" else -1)
        cursor.execute(sql)
        db_connection.commit()            

    cursor.close()

def detect_trees(parameter_set, result_base_path, vhm_input_file, envelope, crs, fint_controller):
    logger = logging.getLogger()

    parameterset_id = parameter_set["id"]
    folder_name = str(parameterset_id)
    output_folder = os.path.join(result_base_path,folder_name)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    print("Detecting " + output_folder ) 
        
    # Configure FintController according to the parameter set
    fint_controller.set_working_dir(output_folder)
    fint_controller.m_altitude_allowed = parameter_set["altitutde_allowed"]
    fint_controller.set_dbh_function(parameter_set["dbh_function"], parameter_set["altitutde_allowed"])
    fint_controller.set_diameter_randomization(parameter_set["randomized"],parameter_set["random_variance"])
    fint_controller.set_minimum_height(parameter_set["minimum_tree_height"])
    fint_controller.set_minimum_detection_height(parameter_set["minimum_detection_tree_height"])
    if parameter_set["output_suffix"]:
        fint_controller.set_output_suffix(parameter_set["output_suffix"])
    else:
        fint_controller.set_output_suffix("")
    if parameter_set["gauss_sigma"]:
        fint_controller.set_gauss_filter(size = parameter_set["gauss_size"], sigma = parameter_set["gauss_sigma"])
    if parameter_set["resize_resolution"]:
        fint_controller.set_resize_resolution(parameter_set["resize_resolution"],parameter_set["resize_method"])

    #Extract VHM
    vhm_output_file = os.path.join(output_folder,"vhm.tif")
    try:
        crop_image(vhm_input_file, vhm_output_file, [envelope]) # Nodata Value may depend on source!
    except ValueError as ve:
        logger.error("Parameterset_ID: "+str(parameterset_id)+" "+result_base_path)
        logger.error(traceback.format_exception(*sys.exc_info()))
        return None

    try:
        #Run FINT
        fint_controller.set_working_dir(output_folder)
        fint_controller.set_normalized_model_file_name(vhm_output_file,None)
        fint_controller.run_process()    
    except Exception as e:
        logger.error("Parameterset_ID: "+str(parameterset_id)+" "+result_base_path)
        logger.error(traceback.format_exception(*sys.exc_info()))
        return None
    finally:
        filelist = [ f for f in os.listdir(output_folder) if f.endswith(".tif") ]
        for f in filelist:
            os.remove(os.path.join(output_folder, f))
        #os.remove(vhm_output_file)

    suffix = fint_controller.m_output_suffix
    treefile_name = "Ind_trees{0}.csv".format("_"+suffix if suffix else "")
    fint_tree_path = os.path.join(output_folder, treefile_name)

    if not os.path.isfile(fint_tree_path):
        return None
        #TODO: Throw Exception?

    df = pd.read_csv(fint_tree_path, delimiter=";", header=None, names=["x","y","hoehe","hoehe_modified","bhd","dominanz"],dtype={"x":np.float32,"y":np.float32,"hoehe":np.str,"hoehe_modified":np.str,"bhd":np.str,"dominanz":np.str})
    geometry = [Point(xy) for xy in zip(df.x, df.y)]
    fint_trees = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    
    return fint_trees 

#Function that returns the 2D euclidean distance  between two trees
def dist_2d(tree1, tree2):
    dx = tree1["x"]-tree2["x"]
    dy = tree1["y"]-tree2["y"]
    d = math.sqrt(dx**2+dy**2)
    return d

def merge_processing(fint_trees_df):
    r_max = 1.5
    dh_max = 5
    #2m oder sigma 2 5/7
    matching_configs = {32: {"Source":32, "Target":31, "rmax":r_max}, 
    #33: {"Source":33, "Target":31, "rmax":r_max}, 
    34: {"Source":34, "Target":31, "rmax":r_max}, 
    #35: {"Source":35, "Target":31, "rmax":r_max}, 
    #36: {"Source":36, "Target":31, "rmax":r_max}
    }

    if len(fint_trees_df)==0:
        return #Nothing to do here


    #prep data
    fint_trees_df["matched32"] = 0
    fint_trees_df["matched33"] = 0
    fint_trees_df["matched34"] = 0
    fint_trees_df["matched35"] = 0
    fint_trees_df["matched36"] = 0

    fint_trees_df["gid"] = fint_trees_df.index

    #ready comparison
    aoi = None #Area of Interest
    dd2d = None
    h_test = None
    dh = None
    dh_max = 5

    #run comparison

    for key, values in fint_trees_df.groupby(["perimeter_id","parameterset_id"]):
        perimeter_id = key[0]
        parameterset_id = key[1]
        if not (parameterset_id in matching_configs):
            continue

        mc = matching_configs[parameterset_id]
        #print("Matching",mc)
        r_max = mc["rmax"]
                
        aoi_trees = fint_trees_df[(fint_trees_df["parameterset_id"]==mc["Target"]) & (fint_trees_df["perimeter_id"]==perimeter_id)]

        if len(aoi_trees) == 0:
            continue 

        aoi_trees_index = aoi_trees.sindex
        
        #Candidate search
        assigned_trees = []
        tree_mappings = {}
        unmapped_trees = []

        sorted_test_trees = values.sort_values(["hoehe"],ascending=False)    
        for tindex, ttree in sorted_test_trees.iterrows():
            #TODO: Check if FINT tree is in Plot bbox? -> Raster was choses a bit larger
            h_test = ttree["hoehe"]
            candidates = []

            x_tt = ttree["x"]
            y_tt = ttree["y"]
            
            #candidate_trees = aoi_trees.cx[x_tt-r_max:x_tt+r_max, y_tt-r_max:y_tt+r_max]        
            candidate_idx = list(aoi_trees_index.intersection([x_tt-r_max,y_tt-r_max,x_tt+r_max,y_tt+r_max]))
            candidate_trees = aoi_trees.iloc[candidate_idx]

            dd2d_min = None
            dh_min = None
            match_candidate = None
            for rindex, rtree in candidate_trees.iterrows():
                #Additionally, already assigned neighboring Reference trees cannot become candidates.
                if rtree.gid in assigned_trees:
                    continue        

                dd2d = dist_2d(ttree,rtree)

                dh = abs(h_test-rtree["hoehe"]) #TODO: Check if abs() is OK

                if dh < dh_max:
                    if (dh_min == None) or (dh<dh_min):
                        dh_min = dh
                        match_candidate = rtree

            if type(match_candidate)!=type(None):
                fint_trees_df.loc[fint_trees_df.gid==match_candidate.gid,"matched"+str(parameterset_id)] = ttree.gid
                assigned_trees.append(match_candidate.gid)                

    #add combined sets to df
    aoi_trees = fint_trees_df[(fint_trees_df["parameterset_id"]==31)]
    aoi_trees_filtered_37 = aoi_trees[(aoi_trees["matched32"]!=0) | (aoi_trees["matched34"]!=0) ].copy()
    aoi_trees_filtered_37["parameterset_id"] = 37

    #aoi_trees_filtered_38 = aoi_trees[(aoi_trees["matched33"]!=0) | (aoi_trees["matched35"]!=0) | (aoi_trees["matched36"]!=0) ].copy()
    #aoi_trees_filtered_38["parameterset_id"] = 38

    #ret = gpd.GeoDataFrame(pd.concat([fint_trees_df,aoi_trees_filtered_37,aoi_trees_filtered_38], ignore_index=True), crs=fint_trees_df.crs)
    ret = gpd.GeoDataFrame(pd.concat([fint_trees_df,aoi_trees_filtered_37], ignore_index=True), crs=fint_trees_df.crs)
    return ret


def process_fst(record,db_connection):
    cursor = db_connection.cursor()



    table_schema = record["table_schema"]
    table_base_name = record["table_base_name"]

    fst_tile_insert_template = "INSERT INTO "+table_schema+"."+table_base_name+"_fst_raster(tile_id, hdom, dg, nh, fst, geom, perimeter_id, flaeche_id) VALUES ({0}, {1}, {2}, {3}, {4}, ST_SetSRID(ST_GeomFromText('{5}'),{6}),{7},{8});"

    result_base_path = record["result_base_path"]
    perimeter_buffer = record["perimeter_buffer"]
    r_max = record["r_max"]
    epsg = record["epsg"]
    crs = record["crs"]

    perimeter_id = record["perimeter_id"]
    flaeche_id = record["flaeche_id"]
    folder_name = "{0}_{1}".format(flaeche_id,perimeter_id)
    output_folder = os.path.join(result_base_path,folder_name)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    print("Processing "+str(perimeter_id))    


    geom = record["geometry"]
    (minx, miny, maxx, maxy) = geom.bounds

    #Envelope by group
    bx = box(minx, miny, maxx, maxy)
    #Envelope by group
    envelope = bx.buffer(perimeter_buffer, resolution=1).envelope

    mischungsgrad = record["mischungsgrad"]
    vhm_input_file_150 = record["vhm_input_file_150"]

    fst_grid_name = "fst_grid.shp"
    fst_grid_path = os.path.join(output_folder,fst_grid_name)

    generate_grid(minx, miny, maxx, maxy, fst_grid_path)

    fst_res = determine_fst(fst_grid_path,vhm_input_file_150,mischungsgrad,envelope)

    if fst_res != 0:
        print("Error calculating FST")
        return
    fst_tiles = gpd.read_file(fst_grid_path)
    fst_tiles["wkt"] = fst_tiles["geometry"].apply(dumps) #reflects onto the underlying df

    for tindex, ttile in fst_tiles[(fst_tiles.intersects(record["geom_flaeche"]))].iterrows():
        (id,hdom,dg,nh,fst) = (ttile["id"],ttile["hdom"],ttile["dg"],ttile["nh"],ttile["FST"])
        if ( 0 <= hdom and hdom <= 100 ):
            sql = fst_tile_insert_template.format(id,hdom,dg,nh,fst,ttile["wkt"],epsg,perimeter_id,flaeche_id)
            cursor.execute(sql)
            db_connection.commit()            

    cursor.close()

def process_trees(record,db_connection):
    cursor = db_connection.cursor()

    logger = logging.getLogger()



    table_schema = record["table_schema"]
    table_base_name = record["table_base_name"]

    fst_select_template = "SELECT gid, geom, fst FROM "+table_schema+"."+table_base_name+"_fst_raster WHERE flaeche_id={0} AND perimeter_id={1}"
    tree_select_template = "SELECT gid, x, y, hoehe, dominanz, bhd, geom, parameterset_id, perimeter_id, flaeche_id, hoehe_modified FROM "+table_schema+"."+table_base_name+"_tree_detected WHERE parameterset_id IN (32,34,37) AND perimeter_id={0} AND flaeche_id={1}"

    tree_insert_template = "INSERT INTO "+table_schema+"."+table_base_name+"_processed_tree(x, y, hoehe, bhd, dominanz, geom, parameterset_id, fst_raster_id, hoehe_modified, fst, flaeche_id) VALUES ({0}, {1}, {2}, {3}, {4}, ST_SetSRID(ST_GeomFromText('{5}'),{6}), {7}, {8},{9},{10},{11});" 
    

    r_max = record["r_max"]
    epsg = record["epsg"]
    crs = record["crs"]
    
    perimeter_id = record["perimeter_id"]
    flaeche_id = record["flaeche_id"]

    print("Processing "+str(perimeter_id))    

    geom = record["geometry"]
    (minx, miny, maxx, maxy) = geom.bounds

    e = len(str(int(minx)))
    f = 10**e
    tile_id = minx*f+miny

    sql_tiles = fst_select_template.format(flaeche_id,perimeter_id)
    fst_tiles_db = gpd.GeoDataFrame.from_postgis(sql_tiles, db_connection, geom_col='geom' )

    if len(fst_tiles_db) == 0:
        return
    

    sql_trees = tree_select_template.format(perimeter_id, flaeche_id)
    fint_trees_db = gpd.GeoDataFrame.from_postgis(sql_trees, db_connection, geom_col='geom' )

    if len(fint_trees_db) == 0:
        return

    fint_trees_db_index = fint_trees_db.sindex
    fst_parameterset_mappings = {111:32,112:32,121:34,122:37,
                                 211:32,212:37,221:34,222:37,
                                 311:37,312:37,321:37,322:37}
    fint_trees_db["fst"] = -1

    picked_trees = []

    for tindex, ttile in fst_tiles_db.iterrows():
        fst = ttile["fst"]
        fst_tile_id = ttile["gid"]
        (tminx, tminy, tmaxx, tmaxy) = ttile["geom"].bounds

        parameterset_id = fst_parameterset_mappings[fst]

        tile_trees = fint_trees_db[(fint_trees_db["parameterset_id"]== parameterset_id) & (fint_trees_db["x"]>=tminx) & (fint_trees_db["x"]<tmaxx) & (fint_trees_db["y"]>=tminy) & (fint_trees_db["y"]<tmaxy)]
        tile_trees["fst_tile_id"] = fst_tile_id
        
        if len(tile_trees)>0:
            tile_trees["fst"] = fst
            picked_trees.append(tile_trees.copy())
        

    if len(picked_trees) == 0:
        return
      
    
    picked_trees_df = gpd.GeoDataFrame(pd.concat(picked_trees, ignore_index=True), crs=fint_trees_db.crs)
    picked_trees_df["geometry"] = picked_trees_df["geom"]
    
    picked_trees_df["wkt"] = picked_trees_df["geometry"].apply(dumps)

    if len(picked_trees_df) == 0:
        return
      
    if "waldmaske_df" in record and type(record["waldmaske_df"])!=type(None):
        waldmaske_df = record["waldmaske_df"]
        pointInForest = sjoin(picked_trees_df, waldmaske_df[[waldmaske_df.geometry.name]], how='left')
        picked_trees_df = pointInForest[np.invert(np.isnan(pointInForest["index_right"]))].drop(columns=["index_right"])

    if perimeter_id == 26331246:
        print("foo")

    if "trasse_mask_df" in record and type(record["trasse_mask_df"])!=type(None):
        trasse_mask_df = record["trasse_mask_df"]
        pointInTrasse = sjoin(picked_trees_df, trasse_mask_df[[trasse_mask_df.geometry.name]], how='left')
        picked_trees_df = pointInTrasse[(np.invert(np.isnan(pointInTrasse["index_right"])) & (pointInTrasse["hoehe"] < 25) ) | (np.isnan(pointInTrasse["index_right"]))]

    for tindex, ttree in picked_trees_df.iterrows():
        sql = tree_insert_template.format(ttree["x"],ttree["y"],ttree["hoehe"],ttree["bhd"],ttree["dominanz"],ttree["wkt"],epsg,ttree["parameterset_id"],ttree["fst_tile_id"],ttree["hoehe_modified"],fst,ttree["flaeche_id"])
        cursor.execute(sql)
        db_connection.commit()   

    cursor.close()




