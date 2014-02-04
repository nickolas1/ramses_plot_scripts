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
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

def _CO13(field, data):
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

infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
filname = 'reduced_'+str(snap).zfill(5)+'/filaments'+str(axis)+'_'+str(snap).zfill(5)+'.txt'

(lmin, lmax) = get_level_min_max(infoname)
(boxlen, unit_l) = get_boxsize(infoname)

ds = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])

add_field('CO13', function=_CO13)

# make a projection of the CO13 density
for i in xrange(1):
    COmap = 'reduced_'+str(snap).zfill(5)+'/CO13sufacedensity_'+str(i)+'.hdf5'
    if not os.path_exists(COmap):
        proj = pf.h.proj('CO13', i)
        frb = proj.to_frb(width, res, center = cntr, height = height)
        f = h5py.File(COmap, 'w')
        dset = f.create_dataset('surface_density', data = np.log10(frb['CO13']))
        f.close()

# read in the rectangles that define the filaments we're interested in
# these are in units of pixels in the finder image, so we will need to translate these
# to unitary units!
rectdata = ascii.read(filname)

for fil in rectdata:
    fildir = 'reduced_'+str(snap).zfill(5)+'/filaments'+str(axis)+'_'+str(fil[0])+'/'
    if not os.path.exists(fildir):
        os.makedirs(fildir)
        
    leftpoint = np.array([fil[1], fil[2]])
    rightpoint = np.array([fil[3], fil[4]])
    width = fil[5]
    # translate to unitary units
    leftpoint /= 2**lmax
    rightpoint /= 2**lmax
    width /= 2**lmax
    print leftpoint
    print rightpoint
    print width
    
    # vector pointing along the filament box's long axis
    vec = rightpoint - leftpoint
    length = np.linalg.norm(vec)
    # the orthogonal direction
    orthovec = (-vec[1], vec[0])
    # normalize them
    orthovec /= np.linalg.norm(orthovec)
    vec /= np.linalg.norm(vec)
    
    # we will move along the lower-left line of the rectangle, starting from startpoint 
    # and ending at endpoint 
    startpoint = leftpoint - orthovec * width/2
    endpoint = rightpoint - orthovec * width/2
    
    # for the distance to move in each step, choose a step size closest to the number of
    # grid cells (at the finest refinement level) it would take to travers the box 
    nl = int(length * 2**lmax)
    nw = int(width * 2**lmax)
    dl = length / nl
    dw = width / nw
    print nl, nw
    print dl, dw
    print length, width, 1/(2**lmax)
    
    print startpoint
    print endpoint
    print dl
    
    vmax = 4.e5
    vmin = -vmax
    vmax = 4.e5
    vmin = -4.e5
    # roughly match hacar et al by takin 0.075 km/s bins
    bins = (vmax - vmin) / 1.e5 / 0.1
    binvals = np.arange(vmin, 1.000001*vmax, (vmax - vmin) / bins)
    binmids = 0.5 * (np.roll(binvals, -1) + binvals)
    binmids = binmids[:len(binmids) - 1]
    # get a version of the bins in km/s instead of cgs
    binmidskms = binmids / 1.e5
    
    f = h5py.File(fildir+'spectrumvels.hdf5', 'w')
    dset = f.create_dataset('binmidskms', data = binmidskms)
    f.close()
    histim = []
    
    for il in xrange(nl):
        print 'step %d of %d along the filament' % (il, nl)
        l = startpoint + vec * dl * (il + 0.5)
        totalhist = np.zeros(bins)
        thesehists = []
        for iw in xrange(nw):
            pt = l + orthovec * dw * (iw + 0.5)
            ray = ds.h.ortho_ray(axis, pt)
            x = ray['x']
            vx = ray['x-velocity']
            dx = ray.fwidth[:,0]
            rho = ray['CO13']
            weight = dx * rho
            #weight /= np.sum(weight)
            """inds = np.argsort(x)
            xs = np.take(x, inds)
            vxs = np.take(vx, inds)
            dxs = np.take(dx, inds)
            rhos = np.take(rho, inds)"""
            hist, binedges = np.histogram(
                vx, 
                range = (vmin, vmax), 
                bins = binvals,
                weights = weight)
            thesehists.append(hist)
            totalhist += hist
        f = h5py.File(fildir+'spectra_'+str(il).zfill(4)+'.hdf5', 'w')
        dset = f.create_dataset('spectra', data = thesehists)
        f.close()
        
        f = h5py.File(fildir+'combinedspectrum_'+str(il).zfill(4)+'.hdf5', 'w')
        dset = f.create_dataset('spectra', data = totalhist)
        f.close()
        
        

    
