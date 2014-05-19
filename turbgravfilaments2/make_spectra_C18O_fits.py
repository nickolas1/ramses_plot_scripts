from __future__ import division


import numpy as np
import sys
import h5py
import shutil
import glob
from astropy.io import fits
from os.path import expanduser

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

snap = int(sys.argv[1])
axis = 2
fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

infofile = fileprefix+'info_'+str(snap).zfill(5)+'.txt'
(time, unit_t) = get_time(infofile)
time_kyr = time * unit_t / (3.15569e7 * 1000)
timestr = str(int(round(time_kyr))).zfill(4)
print timestr


f = h5py.File(fileprefix+'posvel_'+str(axis)+'/spectrumvelsC18O.hdf5')
vels = np.array(f['binmidskms'])
f.close()

outcube = np.zeros([inres, inres, len(vels)], dtype=np.float32)

#
# j is image up
# |
# | 
# |_______i is image right
#
# spectra files march along j
#

### NOTE
# the hdf5 files were created by summing along velocity bins rather than integratig
# along them. we need to multiply them all by 1 / bin width, which is 0.05 km/s!

binwidth = 0.05

filelist = glob.glob(fileprefix+'posvel_'+str(axis)+'/spectra_N2Hplus_*.hdf5')
print filelist[0]

totnonzero = 0
specfile = filelist[0]
f = h5py.File(specfile)
specs = f['spectraC18O']
zerospec = np.zeros(specs.shape)
f.close()

for inj in xrange(inres):
    specfile = fileprefix+'posvel_'+str(axis)+'/spectra_C18O_'+str(inj).zfill(4)+'.hdf5'
    if os.path.isfile(specfile):
        f = h5py.File(specfile)
        specs = np.array(f['spectraC18O']) / binwidth
        print specfile
    else:
        specs = zerospec
    outcube[:, inj, :] = specs
    if os.path.isfile(specfile):
        f.close()
    
hdu = fits.PrimaryHDU(outcube)
hdulist = fits.HDUList([hdu])
hdulist.writeto('spectra.C18O.'+timestr+'.fits')



    

   
