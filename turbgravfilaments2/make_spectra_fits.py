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
inres = 1024

snap = 18
axis = 2
fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

f = h5py.File(fileprefix+'posvel_'+str(axis)+'/spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

outcube = np.zeros([inres, inres, len(vels)])

#
# j is image up
# |
# | 
# |_______i is image right
#
# spectra files march along j
#
totnonzero = 0
for inj in xrange(inres):
    specfile = fileprefix+'posvel_'+str(axis)+'/spectra_C18O_'+str(inj).zfill(4)+'.hdf5'
    f = h5py.File(specfile)
    specs = f['spectraC18O']
    print specfile
    outcube[:, inj, :] = specs
    f.close()
    
hdu = fits.PrimaryHDU(outcube)
hdulist = fits.HDUList([hdu])
hdulist.writeto('spectra_C18O.fits')

    

   
