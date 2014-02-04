from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
from astropy.io import ascii
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font main
    fname=fontdir+'Gotham-Book.ttf', size=13)    
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=11)  


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    
    projaxis = 1

    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)
    # set column density limits so that the images appear the same.
    # low dens cloud has boxsize 10, high dens cloud has boxsize 5.
    # offsets in log of column density limits are thus log10(8)
    if boxlen > 7:
        cdmin = -5.1
        cdmax = -2.0
    else:
        cdmin = -4.19
        cdmax = -1.09

    pf = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    # get a projection of density along y axis
    proj = pf.h.proj('Density', 1)

    wd = 0.625 # this is messed up- figure it out. yt might not get size right.
    # res should be base resolution times 2**levels of refinement * wd
    resx = int(wd * 2**lmax)
    res = (resx,resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    frb = proj.to_frb(width, res, center = cntr, height = height)
    filename = sinkname = 'output_'+str(snap).zfill(5)+'/surface_density_1.hdf5'
    f = h5py.File(filename, 'w')
    dset = f.create_dataset('surface_density', data = np.log10(frb['Density']))
    f.close()
    
    
    # get a projection of density along x axis
    proj = pf.h.proj('Density', 0)
    frb = proj.to_frb(width, res, center = cntr, height = height)
    filename = sinkname = 'output_'+str(snap).zfill(5)+'/surface_density_0.hdf5'
    f = h5py.File(filename, 'w')
    dset = f.create_dataset('surface_density', data = np.log10(frb['Density']))
    f.close()
    
    # get a projection of density along z axis
    proj = pf.h.proj('Density', 2)
    frb = proj.to_frb(width, res, center = cntr, height = height)
    filename = sinkname = 'output_'+str(snap).zfill(5)+'/surface_density_2.hdf5'
    f = h5py.File(filename, 'w')
    dset = f.create_dataset('surface_density', data = np.log10(frb['Density']))
    f.close()    
    
    del(frb)
    del(pf)
    gc.collect()
    
