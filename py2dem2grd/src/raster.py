# File: raster.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
import numpy as np
import gdal
import rasterio
from gdalconst import GA_ReadOnly
#----------------------------------------------------------


#----------------------------------------------------------
# F U N C T I O N    G E T _ R A S T V A L U E    
#----------------------------------------------------------
#
# Finds the raster value for a given row,col in a raster
# result = function(row,col,gdal raster)
#----------------------------------------------------------
def get_rastvalue(col,row,rdata):
    geoTransform = rdata.GetGeoTransform()
    band = rdata.GetRasterBand(1)
    vals = band.ReadAsArray()
    return vals[col][row]
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G E T _ B O U N D I N G B O X
#----------------------------------------------------------
#
# Finds the bounding box (xmin/ymin/xmax/ymax) for a raster
# result = function(gdal raster)
#----------------------------------------------------------
def get_boundingbox(rdata):
    geoTransform = rdata.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * rdata.RasterXSize
    miny = maxy + geoTransform[5] * rdata.RasterYSize
    bbox = [minx, miny, maxx, maxy]
    data = None
    return bbox
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G E T _ N U M R O W C O L              
#----------------------------------------------------------
#
# Finds the total number of rows and columns in a raster
# result = function(gdal raster)
#----------------------------------------------------------
def get_numrowcol(rdata):
    geoTransform = rdata.GetGeoTransform()
    return rdata.RasterXSize, rdata.RasterYSize
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    I S I N R A S T E R              
#----------------------------------------------------------
#
# Finds if a give mesh node is inside a raster
# result = function(x,y,gdal raster)
#----------------------------------------------------------
def isinraster(x,y,rdata):
    bbox = get_boundingbox(rdata)
    if ( (x > bbox[0]) and (x < bbox[2]) and
            (y > bbox[1]) and (y < bbox[3]) ):
        # Inside the raster bbox
        return True
    else:
        return False
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G E T _ R A S T E R S I Z E      
#----------------------------------------------------------
#
# Finds if the cell size of a raster           
# result = function(gdal raster)
#----------------------------------------------------------
def get_rastersize(rdata):
    geoTransform = rdata.GetGeoTransform()
    return geoTransform[1]
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    C O O R D 2 P I X E L            
#----------------------------------------------------------
#
# Finds the row and column for a given x,y pair in a raster
# result = function(x,y,gdal raster)
#----------------------------------------------------------
def coord2pixel(x,y,rdata):
    gt = rdata.GetGeoTransform()
    if (isinraster(x,y,rdata)):
        col = np.int((x - gt[0]) / gt[1])
        row = np.int((gt[3] - y) / -gt[5])
        return col,row
    else:
        return (-1,-1)

#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    P I X E L 2 C O O R D            
#----------------------------------------------------------
#
# Finds the row and column for a given x,y pair in a raster
# result = function(row,col,gdal raster)
#----------------------------------------------------------
def pixel2coord(col,row,rdata):
    xoff, a, b, yoff, d, e = rdata.GetGeoTransform()
    x = a * col + b * row + xoff
    y = d * col + e * row + yoff
    return(x,y)
#----------------------------------------------------------
