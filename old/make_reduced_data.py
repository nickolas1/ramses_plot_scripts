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
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

    if not os.path.exists(fileprefix):
        os.makedirs(fileprefix)
        
    # copy the infofile and sinkfile to the reduced directory 
    shutil.copy(infoname, fileprefix)
    if os.path.exists(sinkname):
        shutil.copy(sinkname, fileprefix)
    
    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)

    pf = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 0.625 # this is messed up- figure it out. yt might not get size right.
    wd = 1.0
    # res should be base resolution times 2**levels of refinement * wd
    resx = int(wd * 2**lmax)
    res = (resx,resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    
 
    for i in range(3):
    # get projection in each direction
        proj = pf.h.proj('Density', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'surface_density_'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Density']))
        f.close()
    # get 5 slices through each direction 
        if i in [1]:
            slc = pf.h.slice(i, 0.01)
            frb = slc.to_frb(width, res, center = cntr, height = height)
            filename = fileprefix+'density_slice_'+str(i)+'_0.hdf5'
            f = h5py.File(filename, 'w')
            dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
            f.close()    
        
        if i in []:
            slc = pf.h.slice(i, 0.25)
            frb = slc.to_frb(width, res, center = cntr, height = height)
            filename = fileprefix+'density_slice_'+str(i)+'_1.hdf5'
            f = h5py.File(filename, 'w')
            dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
            f.close()  
            
        if i in [0, 1, 2]:
            slc = pf.h.slice(i, 0.5)
            frb = slc.to_frb(width, res, center = cntr, height = height)
            filename = fileprefix+'density_slice_'+str(i)+'_2.hdf5'
            f = h5py.File(filename, 'w')
            dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
            f.close()  
        
        if i in [0, 2]:
            slc = pf.h.slice(i, 0.75)
            frb = slc.to_frb(width, res, center = cntr, height = height)
            filename = fileprefix+'density_slice_'+str(i)+'_3.hdf5'
            f = h5py.File(filename, 'w')
            dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
            f.close()  
        
        if i in []:
            slc = pf.h.slice(i, 0.99)
            frb = slc.to_frb(width, res, center = cntr, height = height)
            filename = fileprefix+'density_slice_'+str(i)+'_4.hdf5'
            f = h5py.File(filename, 'w')
            dset = f.create_dataset('volume_density', data = np.log10(frb['Density']))
            f.close()  
        
    if boxlen > 7:
        sphererad = 27.0
    else:
        sphererad = 13.5
    spherevol = 4.0 * np.pi / 3.0 * (sphererad * 3.086e18)**3
    sp = pf.h.sphere(cntr, (sphererad, 'pc'))
        
    nbins = 128
    dmin = 1.e-25
    dmax = 1.e-17
    profilename = fileprefix+'MassAndVolumeInDensityBins.dat'
    profile = BinnedProfile1D(sp,nbins,'Density',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    
    del(frb)
    del(pf)
    del(profile)
    gc.collect()
    
