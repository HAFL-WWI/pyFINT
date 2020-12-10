######################################################################
# This script is based on the software FINT (C++ implementation v1.10 
# from July 2017; (C) ecorisQ - Luuk Dorren, Nicolas Zuanon)
#
# Author: Christoph Schaller, BFH-HAFL, December 2020
#
# Script with type definitions used by FINT.
######################################################################

from enum import Enum
class ModelFileFormatType(Enum):
    ModelFileFormatUndef = 1
    ModelFileFormatAscii = 2
    ModelFileFormatTiff = 3

class MessageType(Enum):
    info = 0
    warning = 1
    error = 2

class SSpatialReference(Enum):
    undefined = 1
    corner = 2
    center = 3

class SpatialReferenceType(Enum):
    spatialReferenceUndefined = 1
    spatialReferenceCorner = 2
    spatialReferenceCenter = 3
    spatialReferenceMaxIndex = 4

class FieldModelDescription:
    nbRows = None
    nbCols = None
    xCoord = None
    yCoord = None
    cellSize = None
    noDataValue = None
    spatialReference = None

    def __init__ (self, rows=0, cols=0, x=0.0, y=0.0, cell_size=-0.0, no_data_value=None,spatial_reference=SpatialReferenceType.spatialReferenceUndefined):
        self.nbRows = rows
        self.nbCols = cols
        self.xCoord = x
        self.yCoord = y
        self.cellSize = cell_size
        self.noDataValue = no_data_value
        self.spatialReference = spatial_reference 

class TreeData:
    m_xCoord = None
    m_yCoord = None
    m_dominance = None
    m_height = None
    m_height_modified = None
    m_altitude = None
    m_diameter = None

    def __init__ (self, x=0.0, y=0.0, d=0, h=-1.0, h_mod=-1.0, a=-0.0):
        self.m_xCoord = x
        self.m_yCoord = y
        self.m_dominance = d
        self.m_height = h
        self.m_height_modified = h_mod
        self.m_altitude = a
        self.m_diameter = -1.0
