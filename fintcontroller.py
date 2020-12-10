######################################################################
# This script is based on the software FINT (C++ implementation v1.10 
# from July 2017; (C) ecorisQ - Luuk Dorren, Nicolas Zuanon)
#
# Author: Christoph Schaller, BFH-HAFL, December 2020
#
# Script with core class of the python FINT implementation.
######################################################################

import os
import sys
import math

import numpy as np

import rasterio
from rasterio import windows

import dominancemask 
import dbhparser 
from typedefs import *

import time
import random

class FintController:

    VERSION_MAJOR = 1
    VERSION_MINOR = 0 #3
    VERSION_YEAR = 2019 #2014
    VERSION_MONTH = "August" #"September"

    HEADER_NCOLS = 0
    HEADER_NROWS = 1
    HEADER_XLL = 2
    HEADER_YLL = 3
    HEADER_CELLSIZE = 4
    HEADER_NODATA = 5
    DEFAULT_NODATAVALUE = -9999.0

    #
    # Private variables
    #

    #Filenames
    m_dem_file_name = None
    m_dsm_file_name = None
    m_nsm_file_name = None
    m_spm_file_name = None

    #Header objects
    m_dem_header = None
    m_dsm_header = None
    m_nsm_header = None
    m_spm_header = None

    #Object repesenting raster source object (rasterio object)
    m_dem_src = None
    m_dsm_src = None
    m_nsm_src = None
    m_spm_src = None

    #Arrays with currently loaded data
    m_dem_data = None
    m_dsm_data = None
    m_nsm_data = None
    m_spm_data = None

    m_dem_nodata_value = None
    m_dsm_nodata_value = None

    #Type of currently used raster (ASCII or TIFF)
    m_model_file_format = None

    #Current working directory for input and output
    m_working_dir = None
    #fieldModelDescription m_nsmHeader #unused

    #minimum Height for considering pisels as tree
    m_minimum_tree_height = None

    #Working variables
    m_max_crown_radius_in_cells = None
    m_mask = None
    m_mu_parser = None
    m_mu_expression = None
    m_mu_height = 0
    m_mu_altitude = 0
    m_diameter_random_range = None

    #Setting variables
    m_altitude_allowed = None
    m_force_file_overriding = None
    m_use_normalized_surface_model_as_input = False
    m_abort_request = False
    m_is_processing = False
    
    #Constructor setting feasible standard values based on setting in C++ client
    def __init__(self):
        self.m_model_file_format = ModelFileFormatType.ModelFileFormatUndef
        self.m_max_crown_radius_in_cells = 0
        self.m_mu_height = 0.
        self.m_mu_altitude = 0.

        self.m_minimum_tree_height = 4
        self.m_mu_parser = dbhparser.DbhParser()
        self.set_dbh_function("H^1.25", False)
        self.m_diameter_random_range = 0 

        self.m_altitude_allowed = False
        self.m_force_file_overriding = True
        self.m_use_normalized_surface_model_as_input = False
        self.m_abort_request = False
        self.m_is_processing = False 

        self.m_mask = dominancemask.DominanceMask()
    
    #Counter Variable and Functions for "progress bar"/progress console output.
    m_max_progress = 0
    def init_progress_bar(self, max_value):
        self.m_max_progress = max_value
        print("Fint v {0}.{1} \n\nStart processing\n".format(self.VERSION_MAJOR, self.VERSION_MINOR))

    def reset_progress_bar(self):
        self.m_max_progress = 0

    def set_progress_bar(self, value): #Print to console instead of graphical progress
        print("\rProcessing {0}%".format(round(value * 100. / self.m_max_progress,2)))

    #Helper hethods for console output
    def display_message (self, message):
        os.system("color 7") #white
        print(message)

    def display_warning (self, warning):
        os.system("color 6") #yellow
        print(warning)

    def display_error (self, error):
        os.system("color 4") #red
        print(error)

    #Set flag for process termination
    def stop_process(self):
        self.m_abort_request = True

    ####
    ## Fintcontrollercore.cpp (mostly)
    ####

    #Check if detection process is running
    def isProcessing(self):
        return self.m_is_processing

    #Set the expression for calculating the DBH based on height.
    def set_dbh_function(self, expression, altitude_allowed ):
        ok = True
        self.m_mu_parser.clear_var()
        self.m_mu_parser.define_var( "H", self.m_mu_height )
    
        self.m_altitude_allowed = altitude_allowed
        if (self.m_altitude_allowed):
            self.m_mu_parser.define_var( "alt", self.m_mu_altitude )
    
        self.m_mu_parser.set_expr(expression)
        try:
          result = self.m_mu_parser.eval()
          self.m_mu_expression = expression
        except Exception as e:
            print(str(e))
            ok = False
        
        return ok
    
    #Run the actual detection process
    def run_process(self):
        self.m_abort_request = False
        self.m_is_processing = True

        ok = self.load_nsm_header()
        if not ok:
            self.terminate_process(1)
            return

        if ( not self.m_use_normalized_surface_model_as_input ):
            ok = self.save_nsm_header()

        if (not ok):
            self.terminate_process(1)
            return

        #Initiate arrays for storing currently processed rows
        self.m_dem_data = []
        self.m_dsm_data = []
        self.m_nsm_data = []
        self.m_spm_data = []

        nRows = self.m_nsm_header.nbRows
        nCols = self.m_nsm_header.nbCols
        cellSize = self.m_nsm_header.cellSize
        trees = []
    
        #// Estimate of the maximum crown diameter (in m; maximum = 30, but for small rasters we have a problem if the recouvrement is bigger than half the raster height or width - 1)
        max_crown_diameter = min( [int( math.floor( min( [nCols, nRows] ) / 2 ) - 1 ), 30] )        
        #// NZ 2015-01-21: need to put a max value of 15 here: if the cellsize is < 1, m_maxCrownRadiusInCells can exceed 15 and therefore the dimensions of the mask
        self.m_max_crown_radius_in_cells = min( [15, int(math.ceil( max_crown_diameter / ( 2 * cellSize ) ) )] )
            

        #// data will be processed block by block. A block has the width of the model and has a height equal to "blockSize"
        #=> the "block" is essentially a rectangular window moving over the raster
        blockSize = self.compute_block_size( nCols )
        if (blockSize==0):
            self.terminate_process(1)
            return

        #Calculate number of blocks to be processed. Mainly used for progress couter
        nbBlocks = int(nRows / blockSize)
        if ( nRows % blockSize > 0 ):
            nbBlocks+=1
    
        currentRow = 0
        rowOffset = 0
        self.init_progress_bar( nbBlocks * 4 + 4 ) #// there are 4 steps per block, plus 4 additional steps at the end
        progress = 0
        #Iterate through blocks
        for iBlock in range(nbBlocks):
            if self.m_abort_request: #Currently unlikely since it is not multithreaded
                break
            loadedRows = 0
            if ( self.m_use_normalized_surface_model_as_input ):
                #Read NSM values from source raster
                loadedRows = self.load_file_data_lines( self.m_nsm_src, self.m_nsm_data, blockSize, nCols, currentRow )
                if ( self.m_altitude_allowed ):
                    loadedDEMRows = self.load_file_data_lines( self.m_dem_src, self.m_dem_data, blockSize, nCols, currentRow )
    
                    assert( loadedDEMRows == loadedRows )

                progress += 1
                self.set_progress_bar( progress )
            else:
                #// Calculate NSM from DEM and DSM
                lastIndex = len(self.m_dem_data)

                #Rad values from DEM and DSM rasters 
                loadedRows = self.load_file_data_lines( self.m_dem_src, self.m_dem_data, blockSize, nCols, currentRow )
                loadedDSMRows = self.load_file_data_lines( self.m_dsm_src, self.m_dsm_data, blockSize, nCols, currentRow )

                progress += 1
                self.set_progress_bar( progress )

                #Calculate NSM calues as DSM-DEM and save the values
                assert( loadedRows == loadedDSMRows )
                if ( loadedRows > 0 ):
                    for row in range(lastIndex,lastIndex+loadedRows,1):
                        self.m_nsm_data.append( np.zeros(self.m_dem_data[0].shape) )
                        for col in range(self.m_nsm_header.nbCols):
                            assert( row < len(self.m_dem_data) and col < len(self.m_dem_data[row]) )
                            dsmValue = self.m_dsm_data[row][col]
                            demValue = self.m_dem_data[row][col]
                            if ( dsmValue == self.m_dsm_nodata_value or demValue == self.m_dem_nodata_value ):
                                self.m_nsm_data[-1][col] = self.m_nsm_header.noDataValue 
                            else:
                                self.m_nsm_data[-1][col] =  dsmValue - demValue 
                    #saveNsmData(m_nsmData.begin() + firstRowToAnalyze, m_nsmData.begin() + lastRowToAnalyze, rowOffset + firstRowToAnalyze );
                    self.save_nsm_data(lastIndex, loadedRows, nCols, currentRow)
    
            #// new raster "species" for future use
            if ( not (self.m_spm_file_name == "" or self.m_spm_file_name == None) ):
                spm_header = FieldModelDescription()
                self.reset_file(self.m_spm_src)
                ok,self.m_spm_src = self.load_file_header( self.m_spm_file_name, spm_header)
                self.load_file_data_lines( self.m_spm_src, self.m_spm_data, blockSize, nCols, currentRow )

            progress += 1
            self.set_progress_bar( progress )
 
            firstRowToAnalyze = 0 if iBlock == 0 else self.m_max_crown_radius_in_cells
            lastRowToAnalyze = len(self.m_nsm_data) if iBlock == nbBlocks - 1 else len(self.m_nsm_data) - self.m_max_crown_radius_in_cells
    
            progress += 1
            self.set_progress_bar( progress )

            #for all rows to analyze
            for row in range(firstRowToAnalyze,lastRowToAnalyze,1):
                #for all columns (0 to raster width)
                for col in range(0,nCols,1):
                    dominance = self.calculDominance( row, col )
                    if ( dominance > 0 ): #Pixel is only a tree if dominance>0
                        assert( not self.m_altitude_allowed or ( row < len(self.m_dem_data) and col < len(self.m_dem_data[row]) ) )
                        altitude = self.m_dem_data[row][col] if self.m_altitude_allowed else 0.0

                        assert( row < len(self.m_nsm_data) and col < len(self.m_nsm_data[row]) )
                        trees.append( TreeData( self.xCoord( col ), self.yCoord( row, rowOffset ), 
                                                dominance,  self.m_nsm_data[row][col], altitude ) )
   
            #// remove the old data in nsm, keep only the minimum data for recouvrement.
            nbOfRowsToRemove = len(self.m_nsm_data) - 2 * self.m_max_crown_radius_in_cells
            if ( nbOfRowsToRemove > 0 ):
                if ( len(self.m_nsm_data) >= nbOfRowsToRemove ):
                    del self.m_nsm_data[:nbOfRowsToRemove]
                if ( len(self.m_dem_data) >= nbOfRowsToRemove ):
                    del self.m_dem_data[:nbOfRowsToRemove]
                if ( len(self.m_dsm_data) >= nbOfRowsToRemove ):
                    del self.m_dsm_data[:nbOfRowsToRemove]
                if ( len(self.m_spm_data) >= nbOfRowsToRemove ):
                    del self.m_spm_data[:nbOfRowsToRemove]
                rowOffset += nbOfRowsToRemove
    
            #// make sure that we update the current row only after loading all files!
            currentRow += loadedRows

            progress += 1
            self.set_progress_bar( progress )
    
        #End iteration through all blocks

        #Calculate diameters an save output files (TXT,CSV and INI)
        if ( not self.m_abort_request ):
            ok = self.compute_all_diameters( trees )
            progress += 1
            self.set_progress_bar( progress )

            self.save_tree_file_txt( trees )
            progress += 1
            self.set_progress_bar( progress )

            self.save_ind_trees_csv( trees )
            progress += 1
            self.set_progress_bar( progress )

            self.save_schema_ini()
            progress += 1
            self.set_progress_bar( progress )
        else:
            self.display_message( "Process aborted by user - no files saved.")

        self.terminate_process(0)

    #Terminate and cleanup the process.
    def terminate_process(self, ret ):
        self.reset_file(self.m_nsm_src)
        self.reset_file(self.m_dem_src)
        self.reset_file(self.m_dsm_src)
        self.reset_file(self.m_spm_src)

        self.m_abort_request = False
        self.m_is_processing = False
        self.reset_progress_bar()
        self.display_message( "Done.")
        return ret
        #?: sys.exit(ret)

    #Get the file format used fpr the models
    def file_format(self):
        return self.m_model_file_format

    #Helper function for computing diameters of detected trees by applying the compiled Expression using the parser. Returns a function reference.
    def compute_diameters(self):
        m_parser = dbhparser.DbhParser()
        m_parser.clear_var()
        m_parser.define_var( "H", self.m_mu_height )
        m_parser.define_var( "alt", self.m_mu_altitude )
    
        m_parser.set_expr(self.m_mu_expression)
        
        def operator(item):
            m_parser.clear_var()
            m_parser.define_var( "H", item.m_height )
            m_parser.define_var( "alt", item.m_altitude )
            item.m_diameter = m_parser.eval()
            return True

        return operator

    #Randomize the passed diameter value within the configured percentage range
    def randomize_diameter (self, percentage):
        m_random = random.Random(time.localtime(0))
        def operator(item):
            item.m_diameter = item.m_diameter * m_random.uniform(1-percentage,1+percentage)
            return True
        return operator

    #Compute the diameter for all detected trees using the configured expression and randomization range
    def compute_all_diameters(self, all_trees):           
        dia_operator = self.compute_diameters()
        for t in all_trees:
            dia_operator(t)

        if ( self.m_diameter_random_range > 0 ):
            rand_operator = self.randomize_diameter(self.m_diameter_random_range / 100.0)
            for t in all_trees:
                rand_operator(t)
        return True

    #Get x coordinate based on passed index x and the raster metadata
    def xCoord(self, x):
        #return self.m_nsm_data.get_tranform[0] + ( x  + ( m_nsmHeader.spatialReference == spatialReferenceCorner ? 0.5 : 0 ) ) * m_nsmHeader.cellSize
        return self.m_nsm_src.xy(0,x)[0]-( 0.5 if self.m_nsm_header.spatialReference == SpatialReferenceType.spatialReferenceCenter else 0) * self.m_nsm_header.cellSize
    
    #Get y coordinate based on passed index y and the raster metadata
    def yCoord(self, y, offset):
        #return self.m_nsm_data.get_tranform[3] + ( self.m_nsm_data.height - y - 1 - offset  + ( self.m_nsm_data.spatialReference == spatialReferenceCorner ? 0.5 : 0 ) ) * self.m_nsm_data.res[0]
        #return self.m_nsm_data.get_tranform[3] + ( self.m_nsm_data.height - y - 1 - offset + self.m_nsm_data.res[0]
        
        #return self.m_nsm_src.xy(self.m_nsm_header.nbRows - y - 1 - offset,0)[1] #TODO: Check validity for FINT logic
        return self.m_nsm_src.xy(y + offset,0)[1]-( 0.5 if self.m_nsm_header.spatialReference == SpatialReferenceType.spatialReferenceCenter else 0) * self.m_nsm_header.cellSize #TODO: Check validity for FINT logic

    #Calculate the Dominance i.e. the core of the detection
    def calculDominance(self, row, col ):
        #// pour chaque sommet de la liste:
        #// * recherche de max local:
        #//   on commence par regarder un carr� 3x3 autour du sommet (puis 4x4... jusqu'� 2*m_maxCrownRadiusInCells x 2*m_maxCrownRadiusInCells)
        #//   on cherche le max en bordure de carr� (i.e., sur un pseudo-cercle) et on compte par ailleurs combien de points du cercle sont < 1
        #//   on en d�duit le ratio nbdecellsbelowmin
        #//   Puis si le max trouv� est < � la hauteur de l'arbre et si le ratio est < � 50%, on accroit la dominance de l'arbre et on agrandit le cercle
        #//   de recherche

        #// Note for future extension: the "species" raster m_spmData can be used here

        assert ( row < len(self.m_nsm_data) and col < len(self.m_nsm_data[row]))
        tree_height = self.m_nsm_data[row][col]
        if ( tree_height < self.m_minimum_tree_height ): # ignore small trees
            return 0

        assert( self.m_max_crown_radius_in_cells < 16 ) ## Can be 30
        dominance = 0
        #Loop through dominance masks with increasing distance
        for i in range(1,self.m_max_crown_radius_in_cells+1,1):
            #Get coordinate pairs for neighbors in current mask
            neighbours_at_distance = self.m_mask.coords(i)

            neighbour_heights = []
            
            #For all neighbors in mask
            for i_neighbour_at_distance in neighbours_at_distance:
                xDistance = i_neighbour_at_distance[0]
                yDistance = i_neighbour_at_distance[1]

                #Translate relative to absolute coordinates
                xIndex = row + xDistance
                yIndex = col + yDistance
                
                if (not ( yIndex < 0 or yIndex >= self.m_nsm_header.nbCols or xIndex < 0 or xIndex >= len(self.m_nsm_data) )):
                    assert( xIndex < len(self.m_nsm_data) and yIndex < len(self.m_nsm_data[xIndex]))

                #Get neighbor height if coordinates are within bounds; 0 otherwise
                neighbour_height = 0 if ( yIndex < 0 
                                          or yIndex >= self.m_nsm_header.nbCols 
                                          or xIndex < 0 
                                          or xIndex >= len(self.m_nsm_data )) \
                                     else self.m_nsm_data[xIndex][yIndex]
                neighbour_heights.append( neighbour_height )          

            highest_neighbour = max(neighbour_heights)
            max_height = highest_neighbour
            number_of_small_neighbours = sum(height < 1 for height in neighbour_heights)

            #// stop if max_height >= tree_height or more than half of the neighbours are < 1
            #// handle special case where 2 neighbours have the same height (interpolation);
            #// in that case, we consider as a tree the cell that is above and/or on the left
            #// TODO: this does not work if more than 2 contiguous cells have the same height!
            ##if ( max_height == tree_height and dominance == 0 and highest_neighbour - neighbour_heights[0] < 4 ):  #TODO: Check logic of last condition
            ##    dominance+=1
            if (max_height == tree_height and dominance == 0):
                dominance+=1
            elif ( max_height >= tree_height or number_of_small_neighbours >= (len(neighbour_heights) / 2) ):
                return dominance
            else:
                dominance+=1
        return dominance

    #Compute the size for the plocks to be read from file.
    #NOTE: The rasterio windowing functions should be able to cope with arbitrary blocksizes. These values from this method have proven to work. 
    # Remains the question, whether the "manual swapping" approach from the C++ version is necessary or whether the rasterio native random access functions are performant enought.   
    def compute_block_size( self, nCols ):  #TODO: Check necessity with respect to rasterio logic
        #// keep the blocksize rather small, so that the UI is frequently updated
        #// from experience, loading ~150k-200k cells at once is acceptable
        block_size = max([self.m_max_crown_radius_in_cells*2+1, int(math.ceil(200000 / nCols))])
        #// special case: if we are dealing with tiled tiffs, the block size is determined by the size of the tiles

        if (self.m_model_file_format == ModelFileFormatType.ModelFileFormatTiff):
            tile_length1 = 0
            tile_length2 = 0
            if (self.m_use_normalized_surface_model_as_input):
                if ( self.m_nsm_src.is_tiled ):
                    tile_length1 = self.m_nsm_src.block_shapes[0][0] 
                if (self.m_altitude_allowed and self.m_dem_src.is_tiled):
                    tile_length2 = self.m_nsm_src.block_shapes[0][0] 
            else:
                if ( self.m_dem_src.is_tiled ):
                    tile_length1 = self.m_dem_src.block_shapes[0][0] 
                if ( self.m_dsm_src.is_tiled ):
                    tile_length2 = self.m_dsm_src.block_shapes[0][0] 

            if ( tile_length1 > 0 and tile_length2 > 0 and tile_length1 != tile_length2 ):
                self.display_error("This version of FINT does not support tiles of different sizes.") 
                block_size = 0
            elif (tile_length1 > 0):
                block_size = tile_length1
            elif (tile_length2 > 0):
                block_size = tile_length2
        
        return block_size


    ####
    ## Fintcontroller.cpp
    ####

    #Set range value for diameter randomization. The value None disables the use of randomization.
    def set_diameter_randomization(self, random, range):
        self.m_diameter_random_range = range if random else -1

    #Set wehether existing files are to be overwritten.
    # NOTE: pyFINT currently ignores this value and alweys overwrites by default. 
    def set_force_file_overriding(self, force):
        self.m_force_file_overriding = force

    #Set whether to use a NSM raster as input or a Combination of DEM and DSM
    def use_normalized_surface_model_as_input(self, use_model):
        self.m_use_normalized_surface_model_as_input = use_model

    #// set the working directory
    #// error out if the folder does not exist
    def set_working_dir(self, working_dir ):
        ok = os.path.isdir(working_dir)
        if ( ok ) :
            self.m_working_dir = working_dir
            os.chdir(self.m_working_dir)
        return ok
    
    #Determine file format based on extension of source file
    def model_file_format_from_file_info(self, filePath ):
        suffix = os.path.splitext(filePath)[1]
        if ( suffix.lower() == ".tif" ):
            return ModelFileFormatType.ModelFileFormatTiff
        elif ( suffix.lower() == ".txt" or suffix.lower() == ".asc" ):
            return ModelFileFormatType.ModelFileFormatAscii
        else:
            return ModelFileFormatType.ModelFileFormatUndef

    #Set the path and type of the DEM and DSM input rasters. Assumes that both rasters have the same type.
    def set_digital_model_file_names(self, dem_file_name, dsm_file_name ):
        ok = os.path.isfile(dem_file_name) and os.path.isfile(dsm_file_name)\
            and self.model_file_format_from_file_info(dem_file_name) == self.model_file_format_from_file_info(dsm_file_name)
        if ( ok ):
            self.m_dem_file_name = dem_file_name
            self.m_dsm_file_name = dsm_file_name
            self.m_model_file_format = self.model_file_format_from_file_info( dem_file_name )
        return ok
    
    #Set the path and type of the NSM input raster
    def set_normalized_model_file_name(self, nsm_file_name, dem_file_name ):
        ok = os.path.isfile(nsm_file_name)
        if ( ok ):
            self.m_nsm_file_name = nsm_file_name
            self.m_model_file_format = self.model_file_format_from_file_info( nsm_file_name )
    
            if ( self.m_altitude_allowed ):
                ok = os.path.isfile(dem_file_name) and\
                   self.model_file_format_from_file_info(dem_file_name) == self.m_model_file_format
                if ( ok ):
                    self.m_dem_file_name = dem_file_name

        return ok

    #Set the path and type of the species input raster
    def set_species_model_file_name(self, spm_file_name ):
        ok = os.path.isfile(spm_file_name) and\
            self.model_file_format_from_file_info(spm_file_name)  == self.m_model_file_format
        if ( ok ):
            self.m_spm_file_name = spm_file_name
        return ok

    #Set the minimum height for pixels to eb considered trees
    def set_minimum_height(self, min_tree_height ):
        self.m_minimum_tree_height = min_tree_height
        return True

    #Initiate loading or input data
    def load_nsm_header(self):
        if ( self.m_use_normalized_surface_model_as_input ):
            return self.load_header_from_normalized_model()
        else:
            return self.load_headers_from_digital_models()

    #Load input NSM raster as well as corresponding metadata.    
    def load_header_from_normalized_model(self):
        ok = False
        self.reset_file(self.m_nsm_src)
        self.m_nsm_header = FieldModelDescription()
        ok,self.m_nsm_src = self.load_file_header( self.m_nsm_file_name, self.m_nsm_header)
    
        if ( ok and  self.m_altitude_allowed ):
            dem_header = FieldModelDescription()
            self.reset_file( self.m_nsm_src)
            self.m_dem_header = FieldModelDescription()
            ok,self.m_dem_src = self.load_file_header( self.m_dem_file_name, dem_header)
            if (ok):
                ok = self.check_headers( dem_header, self.m_nsm_header )
        return ok

    #Load input DEM and DSM rasters as well as corresponding metadata. Set Metadata for NSM accordingly.
    def load_headers_from_digital_models(self):
        ok = True
        dem_header = FieldModelDescription()
        dsm_header = FieldModelDescription()
    
        self.m_nsm_header = FieldModelDescription()

        self.reset_file(self.m_dem_src)
        ok,self.m_dem_src = self.load_file_header( self.m_dem_file_name, dem_header)
    
        if (ok):
            self.reset_file(self.m_dsm_src)
            ok,self.m_dsm_src = self.load_file_header( self.m_dsm_file_name, dsm_header)

        if (ok):
            ok = self.check_headers( dem_header, dsm_header )
    
        if (ok):
            self.m_nsm_header.nbRows = dem_header.nbRows
            self.m_nsm_header.nbCols = dem_header.nbCols
            self.m_nsm_header.xCoord = dem_header.xCoord
            self.m_nsm_header.yCoord = dem_header.yCoord
            self.m_nsm_header.cellSize = dem_header.cellSize
            self.m_nsm_header.spatialReference = dem_header.spatialReference
            self.m_nsm_header.noDataValue = self.m_dem_src.profile["nodata"] if "nodata" in self.m_dem_src.profile else self.DEFAULT_NODATAVALUE
    
            self.m_dem_nodata_value = dem_header.noDataValue
            self.m_dsm_nodata_value = dsm_header.noDataValue

        return ok

    #Open TIFF raster using rasterio and read "header" data from metadata
    def load_file_tiff_header(self, file_name, descr ):
        ok = file_name != "" and file_name != None
        raster_file = None
        
        if ( not ok ):
            self.display_error( "No file name provided!")
    
        if (ok):
            raster_file = rasterio.open(file_name)
            ok = raster_file != None
    
            if (not ok):
                self.display_error( "Failed to load file {0}.".format( file_name ) )

        if (ok):
            imageHeight = raster_file.height
            imageWidth = raster_file.width
    
            descr.nbRows = imageHeight
            descr.nbCols = imageWidth
            descr.spatialReference = SpatialReferenceType.spatialReferenceCenter

            transform = raster_file.get_transform()
            xPixelResolution = transform[1]
            yPixelResolution = -transform[5]

            xCoordUpperLeft = transform[0]
            yCoordUpperLeft = transform[3]

            if (ok):
                ok = self.set_resolution_and_coord_for_tiff( xPixelResolution, yPixelResolution, xCoordUpperLeft, yCoordUpperLeft, imageHeight, descr )
    
        return [ok,raster_file]

    #Set TIFF specific metata in the passed descriptor
    def set_resolution_and_coord_for_tiff(self, xPixelResolution, yPixelResolution, xCoordUpperLeft, yCoordUpperLeft, imageHeight, descr ):
        ok = xPixelResolution == yPixelResolution or xPixelResolution == -yPixelResolution
        if ( not ok ):
            self.display_error( "Incompatible x- and y-resolution in TIFF.")
        else:
            descr.cellSize = xPixelResolution
            descr.xCoord = xCoordUpperLeft
            #descr.yCoord = yCoordUpperLeft - imageHeight * xPixelResolution
            descr.yCoord = yCoordUpperLeft #TODO: Check if this is correct
            descr.noDataValue = self.DEFAULT_NODATAVALUE
        return ok
 
    #Load raster for the goven filename and save metadata in the passed descriptor.
    def load_file_header(self, file_name, descr ):
        if (self.m_model_file_format == ModelFileFormatType.ModelFileFormatAscii):
            return self.load_file_ascii_header( file_name, descr )
        elif (self.m_model_file_format == ModelFileFormatType.ModelFileFormatTiff):
            return self.load_file_tiff_header( file_name, descr )
        else:
            return [False,None]

    #Open ASCII raster using rasterio and read "header" data from metadata
    def load_file_ascii_header(self, file_name, descr ):
        ok = file_name != "" and file_name != None
        raster_file = None

        if ( not ok ):
            self.display_error( "No file name provided!")
    
        if (ok):
            raster_file = rasterio.open(file_name)
            ok = raster_file != None

            if ( not ok ):
                self.display_error("Unable to open file {0} for reading: please check permissions.".format( file_name ) )

        if ( ok ):
            ok = True
            #// create the fieldModelDescription and fill it with the content of the file header
            
            nbCols = raster_file.width
            nbRows = raster_file.height
            noDataValue = raster_file.profile["nodata"] if "nodata" in raster_file.profile else 0.0
            
            
            descr.nbRows = nbRows
            descr.nbCols = nbCols
            descr.spatialReference = SpatialReferenceType.spatialReferenceCenter

            transform = raster_file.get_transform()
            cellSize = transform[1]
            
            xCoord = transform[0]
            yCoord = transform[3]

            spatialRef = SpatialReferenceType.spatialReferenceUndefined
#            if ( keyword.contains( "corner", Qt::CaseInsensitive ) )
#            {
#                spatialRef = spatialReferenceCorner;
#            }
#            else if ( keyword.contains( "center", Qt::CaseInsensitive ) )
#            {
#                spatialRef = spatialReferenceCenter;
#            }

            descr.nbRows = nbRows
            descr.nbCols = nbCols
            descr.xCoord = xCoord
            descr.yCoord = yCoord
            descr.cellSize = cellSize
            descr.noDataValue = noDataValue
            descr.spatialReference = spatialRef

        return [ok,raster_file]


    #Load requested lines from source using rasterio. Append them to the passed data array
    def load_file_data_lines(self, raster_source, data, row_count, col_count, row_index ):
        window = windows.Window(0, row_index, col_count, row_count)
        lines = raster_source.read(1, window=window)
        #if data:
        #    data = np.vstack((data,lines))
        #else:
        #    data = lines
        for l in lines:
            data.append(l)
        return lines.shape[0]

    def reset_file(self, file ):
        if ( file ):
            file.close()
            del file
 
    #Not actually stream based anymore, since rasterio is handling file access
    def reset_stream(self, stream, tiff_stream ):
        if ( stream ):
            stream.close()
            #? del stream
    #        delete *stream;

        if ( tiff_stream ):
            tiff_stream.close()
            #? del tiff_stream
    #        *tiffStream = NULL;



    #####
    ## fintcotrollerchecks.cpp
    #####

#spatialReferenceType fintController::spatialReferenceFromHeader(QStringList &header)
#{
#    spatialReferenceType spatialRef = spatialReferenceUndefined;
#    QString keyword;
#    int xCoord = 0;
#    Q_ASSERT( header.count() == 6 );
#    extractKeywordAndValueFromString( header.at( HEADER_XLL ), keyword, xCoord );
#    Q_ASSERT( keyword.startsWith( "xll", Qt::CaseInsensitive) );
#
#    if ( keyword.contains( "corner", Qt::CaseInsensitive ) )
#    {
#        spatialRef = spatialReferenceCorner;
#    }
#    else if ( keyword.contains( "center", Qt::CaseInsensitive ) )
#    {
#        spatialRef = spatialReferenceCenter;
#    }
#
#    return spatialRef;
#}
#

    #Different sanity checks
    def check_headers(self, model1, model2 ):
        return self.check_grid_sizes(model1,model2) and\
                self.check_spatial_references(model1,model2) and\
                self.check_coordinates(model1,model2) and\
                self.check_cell_sizes(model1,model2)

    def check_grid_sizes(self, model1, model2):
        ok = True
        if ( model1.nbCols != model2.nbCols ):
            self.display_error( "Inconsistent number of columns in data!" )
            ok = False
  
        if ( model1.nbRows != model2.nbRows ):
            self.display_error( "Inconsistent number of rows in data!"  )
            ok = False
        return ok

    def check_spatial_references(self, model1, model2):
        ok = True
        if ( model1.spatialReference != model2.spatialReference ):
            self.display_error( "Elevation and surface models don't use the same spatial reference"  )
            ok = False
        elif ( model1.spatialReference == SpatialReferenceType.spatialReferenceUndefined ):
            self.display_warning("Undefined or unknown spatial reference in file {0}! Using 'corner' as default.".format( self.m_dem_file_name) )
        return ok

    def check_coordinates(self, model1, model2):
        ok = True
        if ( round( model1.xCoord ) != round( model2.xCoord ) or round( model1.yCoord ) != round( model2.yCoord ) ):
            self.display_error( "Elevation and surface models don't use the same coordinates!")
            ok = False
        return ok

    def check_cell_sizes(self, model1, model2):
        ok = True
        if ( model1.cellSize != model2.cellSize ):
            self.display_error("Elevation and surface models don't use the same cellsize!")
            ok = False
        supportedResolution = [ 0.25, 0.5, 1, 1.5, 2 ]
        if not (model1.cellSize in supportedResolution ):
            self.display_error("The cellsize of the input data generates errors.\n"+
                               "Only cellsizes of 0.25, 0.5, 1, 1.5 or 2 m are allowed.")
            ok = False
        return ok

    #// this checking can only be done after loading the data from dem/dsm files!
    def check_data_set_sizes(self):
        assert( len(self.m_dem_data) > 0 and len(self.m_dsm_data) > 0 )
        ok = True
        if ( self.m_nsm_header.nbCols != len(self.m_dem_data[0]) or self.m_nsm_header.nbCols != len(self.m_dsm_data[0]) ):
            self.display_error("Inconsistent number of columns in data!")
            ok = False
        if ( self.m_nsm_header.nbRows != len(self.m_dem_data) or self.m_nsm_header.nbRows != len(self.m_dsm_data) ):
            self.display_error("Inconsistent number of rows in data!")
            ok = False
        return ok

    ######
    ## fintcontrollersave.cpp
    ######

    #Save basic NSM data based on DSM metadata using rasterio 
    def save_nsm_header(self):
        self.reset_file(self.m_nsm_src)
        
        no_data = self.m_dsm_src.nodata        
        out_meta = self.m_dsm_src.meta.copy()
        out_meta.update({"nodata":no_data})

        format_extension = "tif" if self.m_model_file_format == ModelFileFormatType.ModelFileFormatTiff else "asc"

        self.m_nsm_src = rasterio.open(os.path.join(self.m_working_dir, "nsm.{0}".format(format_extension)), "w", **out_meta)
        
        return True

    #Save calculated NSM raster to file using rasterio
    def save_nsm_data(self, firstDataRow, nRows, nCols, firstRow):
        window = windows.Window(0, firstRow, nCols, nRows)

        lines = np.array(self.m_nsm_data[firstDataRow:],np.float32)
        self.m_nsm_src.write(lines, indexes=1, window=window)

        return True


    #Save the detected Trees to a txt. Makes uses up numpy functions.
    def save_tree_file_txt(self, trees ):
        fileName = os.path.join(self.m_working_dir, "treefile.txt")
        
        treeArr = np.array([[tree.m_xCoord, tree.m_yCoord,tree.m_height] for tree in trees])
        if len(treeArr)>0:
            np.savetxt(fileName, treeArr, fmt="%15.2f", delimiter=" ", newline="\n", header="", footer="", comments="# ", encoding=None)
            self.display_message("Saved {0}".format(fileName))
        else:
            open(fileName, 'a').close()
            self.display_message("Saved empty file {0}".format(fileName))

        return True

    #Save the detected Trees to a CSV. Makes uses up numpy functions.
    def save_ind_trees_csv(self, trees ):
        fileName = os.path.join(self.m_working_dir, "Ind_trees.csv")
        
        treeArr = np.array([[tree.m_xCoord, tree.m_yCoord,tree.m_height,tree.m_diameter,tree.m_dominance] for tree in trees])
        if len(treeArr)>0:
            np.savetxt(fileName, treeArr, fmt=["%.2f","%.2f","%.1f","%.1f","%.1i"], delimiter="; ", newline="\n", header="", footer="", comments="# ", encoding=None)
            self.display_message("Saved {0}".format(fileName))
        else:
            open(fileName, 'a').close()
            self.display_message("Saved empty file {0}".format(fileName))

        return True

    #Write INI File. Values are hardcoded
    def save_schema_ini(self):
        fileName = os.path.join(self.m_working_dir, "schema.ini")
        iniFile = open(fileName,"w") 

        iniFile.write("Ind_trees.csv\n")
        iniFile.write("ColNameHeader=False\n")
        iniFile.write("Format=Delimited(;)\n")
        iniFile.write("CharacterSet=OEM)\n")
        iniFile.write("Col1=X Float\n")
        iniFile.write("Col2=Y Float\n")
        iniFile.write("Col3=Treeheight Float\n")
        iniFile.write("Col4=DBH Float\n")
        iniFile.write("Col5=Dominance Float\n")
 
        iniFile.close() 
        self.display_message("Saved {0}".format(fileName))
        return True

