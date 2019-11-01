# File: app.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import PyAdcirc
from pydem2grd.src.interpolate import interpolate 
#----------------------------------------------------------

def run():
    
    #inmeshfile = raw_input('Name of mesh file: ')
    #outmeshfile = raw_input('Name of output mesh file: ')
    #mfac = raw_input('Multiplication factor (e.g. -1): ')
    mfac = -1.0

    #mymesh = PyAdcirc.Mesh(inmeshfile)
    mymesh = PyAdcirc.Mesh('initial-sub.grd')
    print 'Reading mesh...'
    ierr = mymesh.read()
    if ierr == 0:
        exit(ierr)
    print 'Building element table...'
    mymesh.buildElementTable()
  
    print 'Interpolating...'
    intmesh = interpolate(mymesh,'rasterlist.txt',0.5,mfac)
    
    intmesh.write('fort_z.grd')
