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
    mfac = float(raw_input('Multiplication factor (e.g. -1): '))
    '''

    inmeshfile = "NGOM_SACS_Floodplain_v01_flagged.grd"
    outmeshfile = "NGOM_SACS_Floodplain_v01_2xCA_z.grd"
    mfac = float(-1.0)

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
    intmesh = interpolate(mymesh,'rasterlist.txt',0.5,mfac,imethod)
    
    intmesh.write(outmeshfile)
