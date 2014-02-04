from __future__ import division

from yt.mods import *
from yt.config import ytcfg
import gc
import sys
import h5py
import shutil
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

    ds = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 1.0
    # res should be base resolution times 2**levels of refinement * wd
    lmaxplot = min(11, lmax)  
    resx = int(wd * 2**lmaxplot)
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
        del(f)
        del(dset)
        del(frb)
        del(proj)
	gc.collect()
        
    ad = ds.h.all_data()
        
    nbins = 128
    dmin = 1.e-27
    dmax = 1.e-17
    profilename = fileprefix+'MassAndVolumeInDensityBins.dat'
    profile = BinnedProfile1D(ad,nbins,'Density',dmin,dmax,end_collect=True)
    profile.add_fields("CellMassMsun", weight=None)
    profile.add_fields("CellVolume", weight=None)
    profile.write_out(profilename)
    
    del(ds)
    del(ad)
    del(profile)
    gc.collect()
    
