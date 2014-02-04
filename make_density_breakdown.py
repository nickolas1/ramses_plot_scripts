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

hilim = 10**-19.5
lolim = 10**-21.5

hilim = 10**-21.5
lolim = 10**-22.5

def _Highdens(field, data):
    newfield = data['Density']
    antiselection = data['Density'] < hilim
    newfield[antiselection] = 1.e-99
    return newfield

def _Middens(field, data):
    newfield = data['Density']
    antiselection = (data['Density'] < lolim) | (data['Density'] >= hilim)
    newfield[antiselection] = 1.e-99
    return newfield

def _Lowdens(field, data):
    newfield = data['Density']
    antiselection = data['Density'] >= lolim
    newfield[antiselection] = 1.e-99
    return newfield

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
    
    add_field('Highdens', function=_Highdens)
    add_field('Middens', function=_Middens)
    add_field('Lowdens', function=_Lowdens)
    
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
 
    for i in range(1):
    # get projection in each direction
        #proj = pf.h.proj('Density', i)
        proj = pf.h.proj('Highdens', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'surface_density_high_'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Highdens']))
        f.close()
  
        proj = pf.h.proj('Middens', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'surface_density_mid_'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Middens']))
        f.close()
        
        proj = pf.h.proj('Lowdens', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        filename = fileprefix+'surface_density_low_'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Lowdens']))
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
    profilename = fileprefix+'MassAndVolumeInDensityBins_High.dat'
    profile = BinnedProfile1D(sp,nbins,'Highdens',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    del(profile)
    
    profilename = fileprefix+'MassAndVolumeInDensityBins_Mid.dat'
    profile = BinnedProfile1D(sp,nbins,'Middens',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    del(profile)
    
    profilename = fileprefix+'MassAndVolumeInDensityBins_Low.dat'
    profile = BinnedProfile1D(sp,nbins,'Lowdens',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    
    del(frb)
    del(pf)
    del(profile)
    gc.collect()
    
