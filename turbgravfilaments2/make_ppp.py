from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
from os.path import expanduser
from astropy.io import fits


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


# if you were outputting the entire cube, what would the resolution be?
inres = 1024 #64

snap = int(sys.argv[1])
fileprefix = 'output_'+str(snap).zfill(5)+'/'

infofile = fileprefix+'info_'+str(snap).zfill(5)+'.txt'

# left corner of the output box, in units of the cube length
left_corner = np.array([0.5, 0.5, 0.85])
left_corner -= np.mod(left_corner, 1/inres) #align to grid spacing

# output size of the cube in pixels
output_cube_size = 512 #32

outcube = np.zeros([output_cube_size, output_cube_size, output_cube_size], dtype=np.float32)

slicewidth = 1 / inres

ds = load(infofile)

xcenter = left_corner[0] + output_cube_size / 2 / inres
ycenter = left_corner[1] + output_cube_size / 2 / inres

for iz in xrange(output_cube_size):
    zlow = left_corner[2] + iz * slicewidth
    if zlow >= 1:
        zlow -= 1
    zhigh = zlow + slicewidth
    print zlow, zhigh
    
    cr = ds.h.all_data().cut_region(
        ["obj['z'] > "+str(zlow),
         "obj['z'] < "+str(zhigh)])
    
    prj = ds.h.proj('Density', 2, data_source=cr, 
        center = [xcenter, ycenter, (zlow+zhigh)/2])

    frb = prj.to_frb(output_cube_size / inres, output_cube_size) 
    # when storing the projection, get the mean density by dividing by the pixel length
    outcube[:,:,iz] = frb['Density'] / (slicewidth * ds.h.parameter_file.units['cm'])

f = h5py.File('ppp.hdf5', 'w')
dset = f.create_dataset('density', data = outcube, dtype='float32')
f.close()
    
#hdu = fits.PrimaryHDU(outcube)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('spectra.N2H+.'+timestr+'.fits')

    

   
