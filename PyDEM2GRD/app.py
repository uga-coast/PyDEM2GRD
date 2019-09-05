# File: app.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import PyAdcirc
from dem2grd.src.interpolate import interpolate 
#----------------------------------------------------------

def run(): 
    #mymesh = PyAdcirc.Mesh('fort.14')
    mymesh = PyAdcirc.Mesh('fort_refined.grd')
    ierr = mymesh.read()
    if ierr == 0:
        exit(ierr)
    mymesh.buildElementTable()
    #print (mymesh.numElementsAroundNode(1))
    #print (mymesh.elementsAroundNode())
    interpolate(mymesh,'rasterlist.txt')
