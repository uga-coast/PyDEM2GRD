# File: mesh.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
#import pyadcircmodules
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N    G A T H E R V A L U E S      
#----------------------------------------------------------
#
# Sums up the raster pixel values within stencil
# result = function(mesh, raster, values, numvaluesgathered)
#----------------------------------------------------------
def readMesh(meshFileName):
    
    # Open file
    meshFile = open(meshFileName,'r')

    count = 0
    for line in meshFile:
        count += 1

        if count == 1:
            # Read header
            header = line

        elif count == 2:
            tempLine = line.split()
            numElements = int(tempLine[0])
            numNodes = int(tempLine[1])
            print('Number of nodes: ',numNodes)
            print('Number of elements: ',numElements)
    
        elif count < numNodes+3:
            print(line)
        elif count < numNodes+3+numElements:
            print(line)
        else:
            # nodestrings
            print(line)

        #print(line)

    # Read header

    # Close file
    meshFile.close()

    return 1
    #return mesh
