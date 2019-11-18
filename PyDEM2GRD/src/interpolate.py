# File: interpolate.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import PyAdcirc
import rasterio
import sys
import csv
import gdal
import numpy as np
import operator
import math
from shapely.geometry import box    
from shapely.geometry import Point 
from shapely.geometry import Polygon 
from shapely.geometry import mapping
from gdalconst import GA_ReadOnly
from raster import get_rastersize
from raster import get_numrowcol
from raster import get_boundingbox
from raster import pixel2coord
from raster import coord2pixel
#from rasterstats import zonal_stats
from rasterio.mask import mask
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G R I D D A T A
#
#----------------------------------------------------------
def griddata(mesh,meshconn,xc,yc,boundaryNodes,raster,mfac,values,numvaluesgathered):
        
    values = np.zeros(mesh.numNodes())
    numvaluesgathered = np.zeros(mesh.numNodes())
    
    data = gdal.Open(raster, GA_ReadOnly)
    bbox = get_boundingbox(data)
    bboxPoly = box(bbox[0],bbox[1],bbox[2],bbox[3])
  
    # Create a polygon of the Voronoi Diagram about each node
    for i in range(mesh.numNodes()):
        #if mesh.node(i).id() != 14023:
            #continue

        if ( not bboxPoly.contains(Point(mesh.node(i).x(),mesh.node(i).y())) or
                (mesh.node(i).z() > -999.0) ):
            numvaluesgathered[i] = 1
            values[i] = mesh.node(i).z()
            continue

        # Get the number of elements that surround node i
        numElem = mesh.numElementsAroundNode(i)

        # Use the elements that are connected to node i
        # to construct a Voroni Polygon
        pointList = []
        for j in range(numElem):
            # Build a polygon of centroids (vornoi diagram)

            pointList.append((xc[mesh.elementTable(mesh.node(i),j).id()-1],
                yc[mesh.elementTable(mesh.node(i),j).id()-1]))
            
        # Check if the mesh node is a boundary node.
        # If so, then the Voroni Diagram should be constructed with
        # the centroids as well as the mid-point along each boundary
        # line segment and the node coordinates itself
        if (mesh.node(i).id() in boundaryNodes):
            #print 'Boundary node found...'
            # Add current mesh node coordinates to polygon
            pointList.append((mesh.node(i).x(),mesh.node(i).y()))

            if numElem == 1:
                #print mesh.elementTable(mesh.node(i),0).id()
                #print meshconn[mesh.elementTable(mesh.node(i),0).id()-1]
                #n1 = meshconn[mesh.elementTable(mesh.node(i),0).id()-1][0]
                #n2 = meshconn[mesh.elementTable(mesh.node(i),0).id()-1][1]
                #n3 = meshconn[mesh.elementTable(mesh.node(i),0).id()-1][2]
                for k in range(3):
                    if meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k] == mesh.node(i).id():
                        continue

                    neigh = mesh.nodeIndexById(meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k])

                    xmid = 0.5 * (mesh.node(i).x() + mesh.node(neigh).x())
                    ymid = 0.5 * (mesh.node(i).y() + mesh.node(neigh).y())

                    pointList.append((xmid,ymid))

            '''
            # Stil working on this ...
            # For now, just check to see if the current node is in the raster
            else:

                # Find the other line segments that
                # Search the other nodes of the connected elements to find the 
                # other two boundary nodes
                for j in range(numElem):

                    neighEl = mesh.elementTable(mesh.node(i),j).id() # This is an element IDs connected to node i
                    #print neighEl
                    print '' 

                    # Find if any of the nodes that touch the surrounding elements are bondary nodes
                    for k in range(3):
                        if meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k] == mesh.node(i).id():
                            continue
                        #print meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k]
                        
                        if meshconn[neighEl-1][k] == mesh.node(i).id():
                            continue
                        if (meshconn[neighEl-1][k] in boundaryNodes):
                            print meshconn[neighEl-1][k]

                quit()

                    #print meshconn[j] # This is the nodes that make up the element
                    for k in range(3): # 3 b/c 3 nodes make up an element
                        
                        # Skip the current nodes as its already been added to pointList
                        #if meshconn[j][k] == mesh.node(i).id():
                        if meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k] == mesh.node(i).id():
                            continue

                        if (meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k] in boundaryNodes):
                            print 'boundary node'    
                            neigh = mesh.nodeIndexById(meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k])
                            xmid = 0.5 * (mesh.node(i).x() + mesh.node(neigh).x())
                            ymid = 0.5 * (mesh.node(i).y() + mesh.node(neigh).y())
                            
                            pointList.append((xmid,ymid))

                        #####
                        # If one of the connected mesh nodes is on the boundary,
                        # then find the mid-point to the current mesh node
                        # and add that to pointList
                        #if (meshconn[j][k] in boundaryNodes):
                        if (meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k] in boundaryNodes):
                            print meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k]
                            # Find mid-point and add to pointList
                            neigh = mesh.nodeIndexById(meshconn[mesh.elementTable(mesh.node(i),0).id()-1][k])
                            
                            xmid = 0.5 * (mesh.node(i).x() + mesh.node(neigh).x())
                            ymid = 0.5 * (mesh.node(i).y() + mesh.node(neigh).y())
                            #xmid = 0.5 * (mesh.node(i).x() + mesh.node(mesh.nodeIndexById(meshconn[j][k])).x())
                            #ymid = 0.5 * (mesh.node(i).y() + mesh.node(mesh.nodeIndexById(meshconn[j][k])).y())
                            
                            pointList.append((xmid,ymid))
                        ####
        '''

        # Sort the points in a clockwise fashion
        # https://stackoverflow.com/questions/51074984/sorting-according-to-clockwise-point-coordinates
        center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), pointList), [len(pointList)] * 2))
        pointList = (sorted(pointList, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360))
        '''   
        if mesh.node(i).id() == 14023:
            rows = zip(pointList)
            with open('VoroniDiagram.csv', 'wb') as myfile:
                    wr = csv.writer(myfile,sys.stdout, delimiter="\t", quoting = csv.QUOTE_NONE)
                    #wr.writerow(pointList)
                    for row in rows:
                        wr.writerow(row)
            quit()
        ''' 
        
        # Generate a polygon of the pointList
        vor = Polygon(pointList)

        # Mapping converts the vor polygon to a GeoJSON-like mapping of a geometric object.
        # This is needed for the mask operation below
        geoms = [mapping(vor)]
        # https://rasterio.readthedocs.io/en/stable/api/rasterio.mask.html
        # https://rasterio.readthedocs.io/en/stable/topics/masking-by-shapefile.html
        rbbox = []
        with rasterio.open(raster) as src:
            # Build the bounding box of the raster to check if it intersects
            # the stencil for interpolation
            rbbox.append((src.bounds[0],src.bounds[1])) # bottom left
            rbbox.append((src.bounds[0],src.bounds[3])) # top left
            rbbox.append((src.bounds[2],src.bounds[3])) # top right
            rbbox.append((src.bounds[2],src.bounds[1])) # bottom right
            rbbox = Polygon(rbbox)
            # p1.intersects(p2) is true if p1 intersects p2
            if not (rbbox.intersects(vor)):
                continue
            ndv = src.nodata
            out_image, out_transform = mask(src,geoms,crop=True)
        # For Debugging
        # Writes out the masked raster
        # https://rasterio.readthedocs.io/en/stable/topics/masking-by-shapefile.html
        '''
        out_meta = src.meta
        out_meta.update({"driver": "GTiff",
            "height": out_image.shape[1],
            "width" : out_image.shape[2],
            "transform": out_transform})
        with rasterio.open("01_test.tif", "w", **out_meta) as dest:
            dest.write(out_image)
        quit()
        '''

        # Convert the subset matrix to an array
        subset = np.asarray(out_image)
        # Remove no data values to mitigate any overflow issues
        subset = subset[subset > ndv]
        subset = subset * mfac
        numvaluesgathered[i] = subset.size
        values[i] = np.sum(subset)

    return(values,numvaluesgathered)
   
#----------------------------------------------------------
# F U N C T I O N    M E S H C O N N E C T I V I T Y
#----------------------------------------------------------
#
# Finds the connectivitiy of the mesh includes elemental
# connectivity and the boundary nodes
# result = function(mesh)
#----------------------------------------------------------
def meshconnectivity(mesh):

    meshconn = mesh.connectivity()

    # Find the centroid of the elements
    xc = np.zeros(mesh.numElements())
    yc = np.zeros(mesh.numElements())
    for i in range(mesh.numElements()-1):

        x1 = mesh.node(mesh.nodeIndexById(meshconn[i][0])).x()
        x2 = mesh.node(mesh.nodeIndexById(meshconn[i][1])).x()
        x3 = mesh.node(mesh.nodeIndexById(meshconn[i][2])).x()
        
        y1 = mesh.node(mesh.nodeIndexById(meshconn[i][0])).y()
        y2 = mesh.node(mesh.nodeIndexById(meshconn[i][1])).y()
        y3 = mesh.node(mesh.nodeIndexById(meshconn[i][2])).y()

        # Calculate the centroid (xc,yc)
        xc[i] = (x1 + x2 + x3) / 3
        yc[i] = (y1 + y2 + y3) / 3

    # Grab the mesh nodes along the boundary
    boundaryNodesPntr = mesh.boundaryNodes()
    boundaryNodes = []
    for i in range(len(boundaryNodesPntr)):
        boundaryNodes.append(int(boundaryNodesPntr[i].id()))

    return meshconn, xc, yc, boundaryNodes

#----------------------------------------------------------
# F U N C T I O N    G A T H E R V A L U E S      
#----------------------------------------------------------
#
# Sums up the raster pixel values within stencil
# result = function(mesh, raster, values, numvaluesgathered)
#----------------------------------------------------------
def gathervalues(mesh,raster,mfac,values,numvaluesgathered):

    data = gdal.Open(raster, GA_ReadOnly)
    # Find the total number of rows and columns in the raster
    numcols, numrows = get_numrowcol(data)
    # Go ahead and load up raster values
    band = data.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    vals = band.ReadAsArray()
    vals = vals * mfac

    rastersize = get_rastersize(data)
    
    numCells = compute_numcells(mesh,rastersize)
    numCells = np.asarray(numCells)
    N = numCells[0,:]
    CA = numCells[1,:]

    bbox = get_boundingbox(data)
    bboxPoly = box(bbox[0],bbox[1],bbox[2],bbox[3])
    
    for i in range(mesh.numNodes()):

        # Check if node is inside the raster bbox + some buffer
        # Only interpolate on flagged nodes
        #bufr = 2 * N[i] * rastersize
        bufr = 1.2 * N[i] * rastersize
        if ( not bboxPoly.buffer(bufr).contains(Point(mesh.node(i).x(),mesh.node(i).y())) or
                (mesh.node(i).z() > -999.0) ):
            numvaluesgathered[i] = 1
            values[i] = mesh.node(i).z()
            continue

        # Check if the total number of cells have already been acquired
        if (numvaluesgathered[i] == CA[i]):
            continue

        # Check if part of the stencil is inside the raster
        col, row  = coord2pixel(mesh.node(i).x(),mesh.node(i).y(),data)
        
        #print "# Rows, # Columns ",numrows,numcols
        #print "Row, Column ",row, col

        # Get the stencil size in number of rows/cols around the current row/col
        # Left to right and top to bottom nomenclature
        left = int(col - N[i])
        bottom = int(row + N[i])
        right = int(col + N[i])
        top = int(row - N[i])
        
        # Check for negative col/row values in the stencil
        if ( (left < 0) and (bottom < 0) ):
            #print(i+1,'Does not overlap')
            numvaluesgathered[i] = 1
            values[i] = mesh.node(i).z()
            continue

        # Check to make sure the stencil does not go off the raster
        if top < 0:
            top = 0
        if bottom > numrows:
            bottom = numrows
        if left < 0:
            left = 0
        if right > numcols:
            right = numcols

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

            # Sum values in the stencil
            values[i] = values[i] + np.sum(subset)
            numvaluesgathered[i] = numvaluesgathered[i] + np.size(subset)
            '''
            for c in range(left,right):
                for r in range(top,bottom):
                    if ((c > 0) and (c < numcols) and (r < numrows) and (r > 0)):
                        if vals[r][c] > -999:
                            values[i] = values[i] + vals[r][c]
                            numvaluesgathered[i] = numvaluesgathered[i] + 1
            #print 'Total: ',values[i]
            #print 'NumCells: ',numvaluesgathered[i]
            #print ''
            '''
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
def interpolate(mesh,rasterlist,minBathyDepth,mfac,imethod):
    # Grab list of raster files
    f = open(rasterlist,'r')
    files = f.readlines()
    val = np.zeros(mesh.numNodes())
    numval = np.zeros(mesh.numNodes())
    
    meshconn, xc, yc, boundaryNodes = meshconnectivity(mesh)
    
    for f in files:
        # Cycle through each raster
        if imethod == "CA":
            a,b = gathervalues(mesh,f.split()[0],mfac,val,numval)
        elif imethod == "griddata":
            a,b = griddata(mesh,meshconn,xc,yc,boundaryNodes,f.split()[0],mfac,val,numval)
        val = a
        numval = b

    interpvalues = np.zeros(mesh.numNodes())
    for i in range(mesh.numNodes()):
        if (numval[i] == 0):
            interpvalues[i] = mesh.node(i).z()
        elif (numval[i] != 0):
            interpvalues[i] = val[i] / numval[i]
            # Check for minimum bathy depth
        #if interpvalues[i] >= 0 and interpvalues[i] < minBathyDepth:
            #interpvalues[i] = minBathyDepth
    mesh.setZ(interpvalues)
    #mesh.write('fort_z.grd')
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

    # Compute the local mesh size (meters)
    mesh.size = mesh.computeMeshSize()

    sfactor = np.ones(mesh.numNodes())
    for i in range(mesh.numNodes()):
        if (mesh.node(i).z() < -1001 ):
            sfactor[i] = mesh.node(i).z()*-1 - 1000
        else:
            sfactor[i] = 1.0
    
    # Compute N (# of DEM cells radiating omnidirectionally form the cell center)
    N = map(lambda x: (0.25*x)/rastersize, mesh.size)
    N = N * sfactor
    # Compute the total number of DEM cells to average
    N = np.asarray(N)
    N = np.round(N)
    CA = np.piecewise(N, [N < 1, N >= 1], [1, (2*N+1)**2])
    return N, CA
#----------------------------------------------------------
