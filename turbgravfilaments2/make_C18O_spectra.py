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
from scipy import special

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

ds = load(infoname)

# add new density fields
add_field('C18O', function=_C18O)
#add_field('N2Hplus', function=_N2Hplus)

vmax = 2.5e5
vmin = -2.5e5
# roughly match hacar et al by takin 0.05 km/s bins
bins = (vmax - vmin) / 1.e5 / 0.05
binvals = np.arange(vmin, 1.000001*vmax, (vmax - vmin) / bins)
binmids = 0.5 * (np.roll(binvals, -1) + binvals)
binmids = binmids[:len(binmids) - 1]
# get a version of the bins in km/s instead of cgs
binmidskms = binmids / 1.e5
    
# save the velocities to a file
f = h5py.File(specdir+'spectrumvels.hdf5', 'w')
dset = f.create_dataset('binmidskms', data = binmidskms, dtype='float32')
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
    
    rhoC18O = np.array(frb['C18O'])
   # rhoN2Hplus = np.array(frb['N2Hplus'])
    sigmaC18O = 0.0526 # thermal width of C18O line in km/s
    
    sigma = sigmaC18O * 1.e5 # convert to cm/s
    erfdenom = np.sqrt(2*sigma**2)
    
    x = np.array(frb[los])
    dx = np.array(frb[dlos])
    vx = np.array(frb[vlos])
    weight = rhoC18O
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
            # for each point along the slice, march along the projecting dimension
            # and turn each detection into a gaussian. bin this gaussian into the 
            # velbins.
            hist = np.zeros(len(binmids))
            for k in xrange(len(vx[:,i])):
                peak = vx[k,i]
                # calculate the cumulative distribution of this line at each velocity bin edge
                cdfs = 0.5 * (1 + special.erf((binvals - peak) / erfdenom)) * weight[k,i]
                # subtract adjacent values to get the contribution to each bin
                hist = hist + np.diff(cdfs)
            thesehists.append(hist)             

    # once we have the histograms of mass-weighted velocity along each point for this
    # row, save it to an hdf5 file
    f = h5py.File(specdir+'spectra_C18O_'+str(j).zfill(4)+'.hdf5', 'w')
    dset = f.create_dataset('spectraC18O', data = thesehists, dtype='float32')
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
    del(rhoC18O)
    del(weight)
    del(hist)
    del(thesehists)
    gc.collect()
    
    
