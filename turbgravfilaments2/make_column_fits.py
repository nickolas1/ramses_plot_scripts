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
inres = 1024#1024
outres = 256#256
stride = int(inres / outres)
print stride

mu = 2.33

#snap = int(sys.argv[1])
#axis = int(sys.argv[2])

snap = 18
axis = 1

fileprefix = 'reduced_'+str(snap).zfill(5)+'/'


# make the surface density fits files
fileCO = fileprefix+'surface_density_C18O2.hdf5'
fileSD = fileprefix+'surface_density2.hdf5'
fileNH = fileprefix+'surface_density_N2Hplus2.hdf5'

fC18O = h5py.File(fileCO, 'r')
sdC18O = fC18O['surface_density_C18O']
# convert to linear units, divide by mu * mH
sdC18O = 10**np.array(sdC18O) / (mu * const.m_p.cgs.value)
# there are a lot of very small values- cull them before getting some
# info on the range of interesting values  
sdnonzero = sdC18O[sdC18O > 10**4]
print 'mean non-zero C18O column density: ',np.mean(sdnonzero),'cm^-2'
print 'median non-zero C18O column density: ',np.median(sdnonzero),'cm^-2'
print 'max C18O column density: ',np.max(sdC18O),'cm^-2'

hdu = fits.PrimaryHDU(sdC18O)
hdulist = fits.HDUList([hdu])
hdulist.writeto('column_density_C18O.fits')
fC18O.close()


fN2Hplus = h5py.File(fileNH, 'r')
sdN2Hplus = fN2Hplus['surface_density_N2Hplus']
# convert to linear units, divide by mu * mH
sdN2Hplus = 10**np.array(sdN2Hplus) / (mu * const.m_p.cgs.value)
# there are a lot of very small values- cull them before getting some
# info on the range of interesting values  
sdnonzero = sdN2Hplus[sdN2Hplus > 10**4]
print 'mean non-zero N2Hplus column density: ',np.mean(sdnonzero),'cm^-2'
print 'median non-zero N2Hplus column density: ',np.median(sdnonzero),'cm^-2'
print 'max N2Hplus column density: ',np.max(sdN2Hplus),'cm^-2'

hdu = fits.PrimaryHDU(sdN2Hplus)
hdulist = fits.HDUList([hdu])
hdulist.writeto('column_density_N2Hplus.fits')
fN2Hplus.close()


fSD = h5py.File(fileSD, 'r')
sdSD = fSD['surface_density']
# convert to linear units, divide by mu * mH
sdSD = 10**np.array(sdSD) / (mu * const.m_p.cgs.value)
# there are a lot of very small values- cull them before getting some
# info on the range of interesting values  
sdnonzero = sdSD[sdSD > 10**4]
print 'mean non-zero SD column density: ',np.mean(sdnonzero),'cm^-2'
print 'median non-zero SD column density: ',np.median(sdnonzero),'cm^-2'
print 'max SD column density: ',np.max(sdSD),'cm^-2'

hdu = fits.PrimaryHDU(sdSD)
hdulist = fits.HDUList([hdu])
hdulist.writeto('column_density_total.fits')
fSD.close()



   