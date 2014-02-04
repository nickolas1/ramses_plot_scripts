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

def _CO(field, data):
    mu = 2.33 # mean molecular weight
    mH = 1.6733e-24
    lolim = 1000.0 * mu * mH # not interested in anything below 10^3 / cm^3
    hilim = 31622.0 * mu * mH # not interested in anything above 10^4.5 / com^3
    newfield = data['Density']
    antiselection = (data['Density'] < lolim) | (data['Density'] >= hilim)
    newfield[antiselection] = 1.e-99
    return newfield


snap = int(sys.argv[1])
axis = int(sys.argv[2])

if axis == 0:
    los = 'x'
    dlos = 'dx'
    vlos = 'x-velocity'
    sliceax = 'z'
if axis == 1:
    los = 'y'
    dlos = 'dy'
    vlos = 'y-velocity'
    sliceax = 'z'
if axis == 2:
    los = 'z'
    dlos = 'dz'
    vlos = 'z-velocity'
    sliceax = 'y'

infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'

specdir = 'reduced_'+str(snap).zfill(5)+'/posvel_'+str(axis)+'/'
if not os.path.exists(specdir):
    os.makedirs(specdir)

(lmin, lmax) = get_level_min_max(infoname)
(boxlen, unit_l) = get_boxsize(infoname)

ds = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])

add_field('CO', function=_CO)

vmax = 2.5e5
vmin = -2.5e5
# roughly match hacar et al by takin 0.05 km/s bins
bins = (vmax - vmin) / 1.e5 / 0.025
binvals = np.arange(vmin, 1.000001*vmax, (vmax - vmin) / bins)
binmids = 0.5 * (np.roll(binvals, -1) + binvals)
binmids = binmids[:len(binmids) - 1]
# get a version of the bins in km/s instead of cgs
binmidskms = binmids / 1.e5
    
# save the velocities to a file
f = h5py.File(specdir+'spectrumvels.hdf5', 'w')
dset = f.create_dataset('binmidskms', data = binmidskms)
f.close()

"""
    to keep things manageable, make this map on the 1024**3 base grid.
    since refinement is only in regions that are collapsing, and we're
    not interested in those dense regions for the C18O map anyway, this is fine.
"""
res = 2**lmin
dres = 1.0 / res

for j in xrange(200):
    pty = (j + 0.5) * dres
    thesehists = []
    print j, pty
    # get a slice
    slc = ds.h.slice(sliceax, pty)
    
    # get it into a frb
    frb = slc.to_frb(
        (1.0, 'unitary'),           # get the whole extent of the box
        res,                       # don't degrade anything
        center = [0.5, 0.5, 0.5],   # centered in the box
        height = (1.0, 'unitary'))  # get the whole extent of the box
    
    rho = np.array(frb['CO'])
    x = np.array(frb[los])
    dx = np.array(frb[dlos])
    vx = np.array(frb[vlos])
    weight = dx * rho
    # we need to grab rows from the slice differently depending on what axis we're projecting
    if axis == 0:
        for i in xrange(res):
            hist, binedges = np.histogram(
                vx[i,:],
                range = (vmin, vmax),
                bins = binvals,
                weights = weight[i,:])
            thesehists.append(hist)
    if axis > 0:
        for i in xrange(res):
            hist, binedges = np.histogram(
                vx[:,i],
                range = (vmin, vmax),
                bins = binvals,
                weights = weight[:,i])
            thesehists.append(hist) 
                  

    # once we have the histograms of mass-weighted velocity along each point for this
    # row, save it to an hdf5 file
    f = h5py.File(specdir+'spectra_'+str(j).zfill(4)+'.hdf5', 'w')
    dset = f.create_dataset('spectra', data = thesehists)
    dset.attrs['slowindex'] = j
    dset.attrs[sliceax] = pty
    f.close()
    
    del(slc)
    del(frb)
    del(f)
    del(dset)
    del(x)
    del(vx)
    del(dx)
    del(rho)
    del(weight)
    del(hist)
    del(binedges)
    del(thesehists)
    gc.collect()
    
    
