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
    sinkname2 = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.csv'
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

    if not os.path.exists(fileprefix):
        os.makedirs(fileprefix)
        
    # copy the infofile and sinkfile to the reduced directory 
    shutil.copy(infoname, fileprefix)
    if os.path.exists(sinkname):
        shutil.copy(sinkname, fileprefix)
        shutil.copy(sinkname2, fileprefix)
        
    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)

    ds = load(infoname)
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 1.0
    # res should be base resolution times 2**levels of refinement * wd
    resx = int(wd * 2**lmax)
    res = (resx,resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    
 
    for i in range(3):
    # get projection in each direction
        proj = ds.h.proj('Density', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'surface_density_'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Density']))
        f.close()
        del(proj)
        del(frb)
        del(dset)
        gc.collect()
        
    if boxlen > 7:
        sphererad = 27.0
    else:
        sphererad = 13.5
    spherevol = 4.0 * np.pi / 3.0 * (sphererad * 3.086e18)**3
    sp = ds.h.sphere(cntr, (sphererad, 'pc'))
        
    nbins = 128
    dmin = 1.e-25
    dmax = 1.e-17
    profilename = fileprefix+'MassAndVolumeInDensityBins.dat'
    profile = BinnedProfile1D(sp,nbins,'Density',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    
    del(ds)
    del(sp)
    del(profile)
    gc.collect()
    
