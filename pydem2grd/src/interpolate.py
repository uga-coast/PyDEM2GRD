# File: interpolate.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import pyadcircmodules
import rasterio
import sys
import csv
from osgeo import gdal
import numpy as np
import operator
import math
from shapely.geometry import box    
from shapely.geometry import Point 
from shapely.geometry import Polygon 
from shapely.geometry import mapping
from .raster import get_rastersize
from .raster import get_numrowcol
from .raster import get_boundingbox
from .raster import pixel2coord
from .raster import coord2pixel
from rasterio.mask import mask
from functools import reduce
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G A T H E R V A L U E S      
#----------------------------------------------------------
#
# Sums up the raster pixel values within stencil
# result = function(mesh, raster, values, numvaluesgathered)
#----------------------------------------------------------
def gathervalues(mesh,raster,N,CA,mfac,values,numvaluesgathered):

    data = gdal.Open(raster, gdal.GA_ReadOnly)
    # Find the total number of rows and columns in the raster
    numcols, numrows = get_numrowcol(data)
    # Go ahead and load up raster values
    band = data.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    vals = band.ReadAsArray()
    vals = vals * mfac

    rastersize = get_rastersize(data)
    
    bbox = get_boundingbox(data)
    bboxPoly = box(bbox[0],bbox[1],bbox[2],bbox[3])
    
    for i in range(mesh.numNodes()):

        # Check if node is inside the raster bbox + some buffer
        # Only interpolate on flagged nodes
        if mesh.node(i).z() > -999.0:
            continue

        # Check to make sure mesh node is inside the raster
        # plus some buffer = rastersize * N 
        if ( mesh.node(i).x() < bbox[0] - (N[i]+1*rastersize) ) or \
                ( mesh.node(i).x() > bbox[2] + (N[i]+1*rastersize) ) or \
                ( mesh.node(i).y() < bbox[1] - (N[i]+1*rastersize) ) or \
                ( mesh.node(i).y() > bbox[3] + (N[i]+1*rastersize) ):
            continue

        # Check if the total number of cells have already been acquired
        if (numvaluesgathered[i] == CA[i]):
            continue

        # Check if part of the stencil is inside the raster
        col, row, inOut  = coord2pixel(mesh.node(i).x(),mesh.node(i).y(),data)
        
        #print "# Rows, # Columns ",numrows,numcols
        #print "Row, Column ",row, col

        # Get the stencil size in number of rows/cols around the current row/col
        # Left to right and top to bottom nomenclature
        left = int(col - N[i])
        bottom = int(row + N[i])
        right = int(col + N[i])
        top = int(row - N[i])
        
        # Check for negative col/row values in the stencil
        #if ( (left < 0) and (bottom < 0) ):
            #continue
        if top < 0:
            top = 0
        if bottom < 0:
            bottom = 0
        if left < 0:
            left = 0
        if right < 0:
            right = 0

        # Find the x,y coordinates of the stencil and build polygon
        xmin,ymin = pixel2coord(left,bottom,data)
        xmax,ymax = pixel2coord(right,top,data)
        stencilPoly = box(xmin,ymin,xmax,ymax)
        
        if (not stencilPoly.touches(bboxPoly)):
            # Find which stencil cells are contained in the raster
            subset = vals[top:bottom,left:right]
            
            # Convert the subset matrix to an array
            subset = np.asarray(subset)

            # Remove no data values to mitigate any overflow issues
            subset = subset[subset >= -999]
            subset = subset[subset <= 999]
            
            # Check to make sure there are valid elevations in the stencil
            if subset.size == 0:
                continue

            # Check for vertical/raised feature nodes
            if mesh.node(i).z() == -2000:

                # Find the mean and standard deviation
                mean = np.average(subset)
                std = np.std(subset)

                # If the node is in the bounds of more than 1 raster,
                # then keep the highest (minimum for ADCIRC) value.
                if ( (mean - 2*std) < np.min(subset) ):

                    values[i] = np.min(subset)
                    numvaluesgathered[i] = -2000

                else:

                    subset = subset[subset <= (mean - 2*std)]
                    values[i] = np.mean(subset)
                    numvaluesgathered[i] = -2000
                
                print(('Mesh Node ',mesh.node(i).id(),' is a raised feature node w/ elevation: ',values[i]))
                
            # Not flagged as a vertical/raised feature node
            else:

                # Sum values in the stencil
                values[i] = values[i] + np.sum(subset)
                numvaluesgathered[i] = numvaluesgathered[i] + np.size(subset)
        
        else:
            
            print((i+1,'Does not overlap'))
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
def interpolate(mesh,rasterlist,minBathyDepth,mfac,imethod):
    # Grab list of raster files
    f = open(rasterlist,'r')
    files = f.readlines()
    val = np.zeros(mesh.numNodes())
    numval = np.zeros(mesh.numNodes())
   
    # Compute the local mesh size (meters)
    mesh.size = mesh.computeMeshSize()

    rastersize = 0
    for f in files:
        print(f)
        # Cycle through each raster
        # Check to see if raster size changed
        data = gdal.Open(f.split()[0], gdal.GA_ReadOnly)
        newrastersize = get_rastersize(data)
        if (abs(rastersize-newrastersize) > 0.10): # Raster size changed > 10 cm
            # Re-calculate N and CA based on the updated raster size
            print(('Raster size changed from ',rastersize,' to ',newrastersize,'. Re-calculating N & CA.'))
            rastersize = newrastersize
            numCells = compute_numcells(mesh,rastersize)
            numCells = np.asarray(numCells)
            N = numCells[0,:]
            CA = numCells[1,:]

        a,b = gathervalues(mesh,f.split()[0],N,CA,mfac,val,numval)
        
        val = a
        numval = b

    interpvalues = np.zeros(mesh.numNodes())
    for i in range(mesh.numNodes()):

        if (numval[i] == 0):

            interpvalues[i] = mesh.node(i).z()

        elif (numval[i] == -2000):

            interpvalues[i] = val[i]

        elif (numval[i] != 0):

            interpvalues[i] = val[i] / numval[i]

        # Check for minimum bathy depth
        #if interpvalues[i] >= 0 and interpvalues[i] < minBathyDepth:
            #interpvalues[i] = minBathyDepth
    
    mesh.setZ(interpvalues)
    
    return mesh
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
    #mesh.reproject(26916)

    sfactor = np.ones(mesh.numNodes())
    for i in range(mesh.numNodes()):
        
        if (mesh.node(i).z() < -1001 ) and (mesh.node(i).z() > -1100) :
            sfactor[i] = mesh.node(i).z()*-1 - 1000
        
        # values of -2000 or less are flagged as vertical/raised feature nodes
        elif (mesh.node(i).z() == -2000):
            sfactor[i] = 2.0
        
        else:
            sfactor[i] = 1.0
    
    # Compute N (# of DEM cells radiating omnidirectionally form the cell center)
    N = [(0.25*x)/rastersize for x in mesh.size]
    N = N * sfactor
    # Compute the total number of DEM cells to average
    N = np.asarray(N)
    N = np.round(N)
    CA = np.piecewise(N, [N < 1, N >= 1], [1, (2*N+1)**2])
    return N, CA
#----------------------------------------------------------
