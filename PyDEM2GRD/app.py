# File: app.py
# Name: Matthew V. Bilskie

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import PyAdcirc
from pydem2grd.src.interpolate import interpolate 
from pydem2grd.src.interpolate import griddata
#----------------------------------------------------------

def run():
    ''' 
    inmeshfile = raw_input('Name of mesh file: ')
    outmeshfile = raw_input('Name of output mesh file: ')
    rlistfile = raw_input('Name of raster list file: ')
    mfac = float(raw_input('Multiplication factor (e.g. -1): '))
    '''

    inmeshfile = "holes_working_z.grd"     
    outmeshfile = "holes_working_z2.grd"
    mfac = float(-1.0)
    #minbath = 0.25
    #minbath = -0.25
    minbath = 0.00
    rlistfile = 'rasterlist_NED.txt'

    mymesh = PyAdcirc.Mesh(inmeshfile)
    print 'Reading mesh...'
    ierr = mymesh.read()
    if ierr == 0:
        exit(ierr)
    print 'Building element table...'
    mymesh.buildElementTable()
   
    print 'Interpolating...'
    #imethod = 'griddata'
    imethod = 'CA'
    intmesh = interpolate(mymesh,rlistfile,minbath,mfac,imethod)
    
    print  'Saving mesh file...'
    intmesh.write(outmeshfile)

    print 'Finished! :)'
