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

    inmeshfile = "NGOM_SACS_Floodplain_v02_cleaned_flagged.grd"
    outmeshfile = "NGOM_SACS_Floodplain_v02_cleaned_2xCA_z.grd"
    mfac = float(-1.0)
    rlistfile = 'rasterlist_all.txt'

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
    intmesh = interpolate(mymesh,rlistfile,0.5,mfac,imethod)
    
    print  'Saving mesh file...'
    intmesh.write(outmeshfile)

    print 'Finished! :)'
