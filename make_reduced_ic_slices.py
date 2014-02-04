from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
import shutil
from astropy.io import ascii
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    
    fileprefix = 'ic_slices_reduced_'+str(snap).zfill(5)+'/'

    if not os.path.exists(fileprefix):
        os.makedirs(fileprefix)
        
    # copy the infofile and sinkfile to the reduced directory 
    shutil.copy(infoname, fileprefix)
    
    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)

    pf = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 0.625 # this is messed up- figure it out. yt might not get size right.
    wd = 0.5
    # res should be base resolution times 2**levels of refinement * wd
    resx = int(wd * 2**lmax)
    res = (resx,resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    
    i = 0 # get slices along x axis
    dslice = wd/resx
    counter = 0
    for slice in np.arange(cntr[i] - wd/2, cntr[0] + wd/2 + dslice, dslice):
        print counter, slice
        counter += 1
        slc = pf.h.slice(i, slice)
        frb = slc.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'density_slice_'+str(i)+'_'+str(counter).zfill(5)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
        f.close()    
    
        del(slc)
        del(frb)
        gc.collect()
        if counter == 7: sys.exit()
    
