# File: interpolate.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import PyAdcirc
import gdal
import numpy as np
from shapely.geometry import box    
from shapely.geometry import Point 
from gdalconst import GA_ReadOnly
from raster import get_rastersize
from raster import get_numrowcol
from raster import get_boundingbox
from raster import pixel2coord
from raster import coord2pixel
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G A T H E R V A L U E S      
#----------------------------------------------------------
#
# Sums up the raster pixel values within stencil
# result = function(mesh, raster, values, numvaluesgathered)
#----------------------------------------------------------
def gathervalues(mesh,raster,values,numvaluesgathered):

    data = gdal.Open(raster, GA_ReadOnly)
    # Find the total number of rows and columns in the raster
    numcols, numrows = get_numrowcol(data)
    # Go ahead and load up raster values
    band = data.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    vals = band.ReadAsArray()

    rastersize = get_rastersize(data)
    
    numCells = compute_numcells(mesh,rastersize)
    numCells = np.asarray(numCells)
    N = numCells[0,:]
    CA = numCells[1,:]

    bbox = get_boundingbox(data)
    bboxPoly = box(bbox[0],bbox[1],bbox[2],bbox[3])
    
    for i in range(mesh.numNodes()):
        # Check if the total number of cells have already been acquired
        if (numvaluesgathered[i] == CA[i]):
            continue

        # Check if node is inside the raster bbox + some buffer
        bufr = 2 * N[i] * rastersize
        if (not bboxPoly.buffer(bufr).contains(Point(mesh.node(i).x(),mesh.node(i).y()))):
            continue

        # Check if part of the stencil is inside the raster
        col, row  = coord2pixel(mesh.node(i).x(),mesh.node(i).y(),data)

        # Get the stencil size in number of rows/cols around the current row/col
        # Left to right and top to bottom nomenclature
        left = int(col - N[i])
        bottom = int(row + N[i])
        right = int(col + N[i])
        top = int(row - N[i])
        # Check for negative col/row values in the stencil
        if ( (left < 0) and (top < 0) ):
            #print(i+1,'Does not overlap')
            continue

        # Find the x,y coordinates of the stencil and build polygon
        xmin,ymin = pixel2coord(left,bottom,data)
        xmax,ymax = pixel2coord(right,top,data)
        stencilPoly = box(xmin,ymin,xmax,ymax)

        if (not stencilPoly.touches(bboxPoly)):
            #print(i+1,'Overlaps')
            # Find which stencil cells are contained in the raster
            for c in range(left,right):
                for r in range(top,bottom):
                    if ((c > 0) and (c < numcols) and (r < numrows) and (r > 0)):
                        if vals[r][c] > -9999:
                            values[i] = values[i] + vals[r][c]
                            numvaluesgathered[i] = numvaluesgathered[i] + 1
        else:
            print(i+1,'Does not overlap')
            continue

    return(values,numvaluesgathered)
#----------------------------------------------------------
    

#----------------------------------------------------------
# F U N C T I O N    I N T E R P O L A T E        
#----------------------------------------------------------
#
# Cycle through a list of rasters and interpolate
# DEM values to the mesh.
# result = function(mesh, rasterlist)
#----------------------------------------------------------
def interpolate(mesh,rasterlist):
    # Grab list of raster files
    f = open(rasterlist,'r')
    files = f.readlines()
    val = np.zeros(mesh.numNodes())
    numval = np.zeros(mesh.numNodes())
    for f in files:
        # Cycle through each raster
        a,b = gathervalues(mesh,f.split()[0],val,numval)
        val = a
        numval = b

    interpvalues = np.zeros(mesh.numNodes())
    for i in range(mesh.numNodes()):
        #print(numval[i])
        if (numval[i] != 0):
            interpvalues[i] = val[i] / numval[i]
    mesh.setZ(interpvalues)
    mesh.write('fort_z.grd')

    return
#----------------------------------------------------------


#----------------------------------------------------------
# F U N C T I O N    C O M P U T E _ N U M C E L L S
#----------------------------------------------------------
#
# Compute the total numhber of DEM cells that should be 
# interpolated for each mesh node using the CCA method
# of Bilskie & Hagen (2013)
# result = function(mesh, rastersize)
#----------------------------------------------------------
def compute_numcells(mesh,rastersize):

    # mesh -> mesh object
    # rastersize -> floating point of DEM cell size

    # Reproject the mesh to UTM coordinates
    mesh.reproject(26916)

    # Compute the local mesh size (meters)
    mesh.size = mesh.computeMeshSize()

    # Compute N (# of DEM cells radiating omnidirectionally form the cell center)
    N = map(lambda x: (0.25*x)/rastersize, mesh.size)
    # Compute the total number of DEM cells to average
    N = np.asarray(N)
    N = np.round(N)
    CA = np.piecewise(N, [N < 1, N >= 1], [1, (2*N+1)**2])
    return N, CA
#----------------------------------------------------------
