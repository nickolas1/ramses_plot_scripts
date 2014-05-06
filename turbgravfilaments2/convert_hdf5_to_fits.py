from __future__ import division


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
import shutil
import glob
from astropy.io import ascii
from astropy.io import fits
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from scipy import special
from mpl_toolkits.mplot3d import Axes3D

"""
usage:
python fooppv.py N A F
N: number of output to use. reduced_N needs to be here.
A: axis of the projection (0, 1, 2)
"""

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *


# downsample the spectra to a 256 squared grid
inres = 1024 #1024
outres = 1024 #256
stride = int(inres / outres)
print stride


snapstart = int(sys.argv[1])
snapend = int(sys.argv[2])

for snap in xrange(snapstart, snapend, 1):
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    infofile = fileprefix+'info_'+str(snap).zfill(5)+'.txt'
    (time, unit_t) = get_time(infofile)
    time_kyr = time * unit_t / (3.15569e7 * 1000)

    timestr = str(int(round(time_kyr))).zfill(4)
    print timestr
    
    for axis in xrange(3):
        ax = str(axis)   
        # make the surface density fits files
        fileCO = fileprefix+'surface_density_C18O'+ax+'.hdf5'
        fileSD = fileprefix+'surface_density'+ax+'.hdf5'
        fileNH = fileprefix+'surface_density_N2Hplus'+ax+'.hdf5'
        
        outCO = 'columndensity.total.'+timestr+'.'+ax+'.fits'
        outSD = 'columndensity.C18O.'+timestr+'.'+ax+'.fits'
        outNH = 'columndensity.N2H+.'+timestr+'.'+ax+'.fits'
        
        print outCO
        print outSD
        print outNH
        
        fC18O = h5py.File(fileCO, 'r')
        sdC18O = fC18O['surface_density_C18O']
        # convert to linear units
        sdC18O = 10**np.array(sdC18O)
        print 'max C18O column density: ',np.max(sdC18O),'cm^-2'
        
        hdu = fits.PrimaryHDU(sdC18O)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(outCO)
        fC18O.close()


        fN2Hplus = h5py.File(fileNH, 'r')
        sdN2Hplus = fN2Hplus['surface_density_N2Hplus']
        # convert to linear units
        sdN2Hplus = 10**np.array(sdN2Hplus)
        print 'max N2Hplus column density: ',np.max(sdN2Hplus),'cm^-2'

        hdu = fits.PrimaryHDU(sdN2Hplus)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(outNH)
        fN2Hplus.close()


        fSD = h5py.File(fileSD, 'r')
        sdSD = fSD['surface_density']
        # convert to linear units
        sdSD = 10**np.array(sdSD)
        sdnonzero = sdSD[sdSD > 10**4]
        print 'mean non-zero SD column density: ',np.mean(sdnonzero),'cm^-2'
        print 'median non-zero SD column density: ',np.median(sdnonzero),'cm^-2'
        print 'max SD column density: ',np.max(sdSD),'cm^-2'

        hdu = fits.PrimaryHDU(sdSD)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(outSD)
        fSD.close()
        


   