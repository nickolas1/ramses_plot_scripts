from __future__ import division

import gc
import sys
import h5py
from os.path import expanduser
from astropy.io import fits


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *


snap = int(sys.argv[1])

infofile = 'reduced_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'

(time, unit_t) = get_time(infofile)
time_kyr = time * unit_t / (3.15569e7 * 1000)

timestr = str(int(round(time_kyr))).zfill(4)


f = h5py.File('ppp'+str(snap).zfill(5)+'.hdf5', 'r')
outcube = f['density']

hdu = fits.PrimaryHDU(outcube)
hdulist = fits.HDUList([hdu])
hdulist.writeto('ppp.'+timestr+'.fits')

f.close()