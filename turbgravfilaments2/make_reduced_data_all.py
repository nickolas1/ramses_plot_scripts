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

def _C18O(field, data):
    mu = 2.33 # mean molecular weight
    mH = 1.6733e-24
    lolim = 1000.0 * mu * mH # not interested in anything below 10^3 / cm^3
    hilim = 31622.0 * mu * mH # not interested in anything above 10^4.5 / com^3
    newfield = data['Density']
    antiselection = (data['Density'] < lolim) | (data['Density'] >= hilim)
    newfield[antiselection] = 0.0
    return newfield
    
def _N2Hplus(field, data):
    mu = 2.33 # mean molecular weight
    mH = 1.6733e-24
    lolim = 31622.0 * mu * mH # not interested in anything below 10^4.5 / cm^3
    newfield = data['Density']
    antiselection = (data['Density'] < lolim)
    newfield[antiselection] = 0.0
    return newfield    

MakeDensityPDF = False

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    sinkname2 = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.csv'
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

    # copy info files to the reduced directory
    if not os.path.exists(fileprefix):
        os.makedirs(fileprefix)       
        # copy the infofile and sinkfile to the reduced directory 
        shutil.copy(infoname, fileprefix)
        if os.path.exists(sinkname):
            shutil.copy(sinkname, fileprefix)
        if os.path.exists(sinkname2):
            shutil.copy(sinkname2, fileprefix)    
    
    # figure out resolution and box size
    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)

    ds = load(infoname)
    
    # add new density fields
    add_field('C18O', function=_C18O)
    add_field('N2Hplus', function=_N2Hplus)
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 1.0
    # full res should be base resolution times 2**levels of refinement * wd
    # we're going to reduce to a 1024 squared image
    lmaxplot = min(10, lmax)  
    resx = int(wd * 2**lmaxplot)
    res = (resx, resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')

    for i in range(3):
    # get projections in each direction
        proj = ds.h.proj('C18O', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        print 'done creating C18O frb ',i    
        filename = fileprefix+'surface_density_C18O'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density_C18O', data = np.log10(frb['C18O']), dtype='float32')
        print 'done creating HDF5 dset for C18O ',i
        f.close()
        del(proj)
        del(frb)
        del(f)
        del(dset)
        gc.collect() 
         
        proj = ds.h.proj('N2Hplus', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        print 'done creating N2Hplus frb ',i    
        filename = fileprefix+'surface_density_N2Hplus'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density_N2Hplus', data = np.log10(frb['N2Hplus']), dtype='float32')
        print 'done creating HDF5 dset for N2Hplus ',i
        f.close()
        del(proj)
        del(frb)
        del(f)
        del(dset)
        gc.collect()
                
        proj = ds.h.proj('Density', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        print 'done creating Density frb ',i    
        filename = fileprefix+'surface_density'+str(i)+'.hdf5'
        f = h5py.File(filename, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['Density']), dtype='float32')
        print 'done creating HDF5 dset for Density ',i
        f.close()
        del(proj)
        del(frb)
        del(f)
        del(dset)
        gc.collect()  
        
       
    if MakeDensityPDF: 
        ad = ds.h.all_data()
        nbins = 128
        dmin = 1.e-27
        dmax = 1.e-17
        profilename = fileprefix+'MassAndVolumeInDensityBins.dat'
        profile = BinnedProfile1D(ad,nbins,'Density',dmin,dmax,end_collect=True)
        profile.add_fields("CellMassMsun", weight=None)
        profile.add_fields("CellVolume", weight=None)
        profile.write_out(profilename)
        del(ad)
        del(profile)
                       
    del(ds)
    gc.collect()
    
