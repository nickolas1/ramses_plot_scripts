from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
import copy
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
#add_field('C18O', function=_C18O)
add_field('N2Hplus', function=_N2Hplus)

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
    to keep things manageable, make this output map on the 1024**3 base grid.
    outres: the output resolution.
    inres: the resolution we sample the grid at.
"""
outres = 2**lmin
outdres = 1.0 / outres

refinefac = 8
inres = outres * refinefac
indres = 1.0 / inres

# tabulate the error function from -3 to 3 
#erfx = np.arange(-3,3.0001,10/512.)
#erfy = special.erf(erfx)
#def closest_erf_value(xvals, yvals, inval):
#    return yvals[(np.abs(xvals-inval)).argmin()]

fileN2Hplus = 'reduced_'+str(snap).zfill(5)+'/surface_density_N2Hplus'+str(axis)+'.hdf5'
print fileN2Hplus
fN2Hplus = h5py.File(fileN2Hplus, 'r')
sdN2Hplus = fN2Hplus['surface_density_N2Hplus']
# convert to linear units, divide by mu * mH
sdN2Hplus = 10**np.array(sdN2Hplus) / (2.33 * 1.66e-24)

for j in xrange(256):
    print 'n2h+ density sum: ',np.min(sdN2Hplus[j,:]),sdN2Hplus[j,:].sum(),np.max(sdN2Hplus[j,:])
    #for i in xrange(1024): print i,sdN2Hplus[j,i]
    if(sdN2Hplus[j,:].sum() == 0):
        print 'skipping ',j
        continue
    outpty = (j + 0.5) * outdres
    thesehists = []
    print j, outpty
    
    jj = 0
    for rj in xrange(refinefac):
        inpty = (j*refinefac + jj + 0.5) * indres
        print 'inpty: ',inpty
        # get a slice
        slc = ds.h.slice(sliceax, inpty)
    
        # get it into a frb
        frb = slc.to_frb(
            (1.0, 'unitary'),           # get the whole extent of the box
            inres,                       # don't degrade anything
            center = [0.5, 0.5, 0.5],   # centered in the box
            height = (1.0, 'unitary'))  # get the whole extent of the box
    
       # rhoC18O = np.array(frb['C18O'])
        rhoN2Hplus = np.array(frb['N2Hplus'])
       # sigmaC18O = 0.0526 # thermal width of C18O line in km/s
        sigmaN2Hplus = 0.0535 # thermal width of N2Hplus line in km/s
    
       # sigma = sigmaC18O * 1.e5 # convert to cm/s
        sigma = sigmaN2Hplus * 1.e5 # convert to cm/s
        erfdenom = np.sqrt(2*sigma**2)
    
        x = np.array(frb[los])
        vx = np.array(frb[vlos])
        dx = np.array(frb[dlos])
        mindx = np.min(dx)
        print 'max(dx), min(dx), outdres: ',np.max(dx),np.min(dx),outdres
        # weight = dx * rhoC18O
        weight = indres * rhoN2Hplus  # the dx * rho line above isn't necessary for this frb scheme- all cells are the same size, so set to indres
        # we need to grab rows from the slice differently depending on what axis we're projecting
        hist = np.zeros(len(binmids))
        erfvals = np.zeros(len(binvals))
        if axis == 0:
            for i in xrange(inres):
                hist, binedges = np.histogram(
                    vx[i,:],
                    range = (vmin, vmax),
                    bins = binvals,
                    weights = weight[i,:])
                thesehists.append(hist)
        if axis > 0:
            i = 0
            for ii in xrange(inres):
                # for each point along the slice, march along the projecting dimension
                # and turn each detection into a gaussian. bin this gaussian into the 
                # velbins.
                hist[:] = 0
                k = 0
                dkmin = 1
                for ik in xrange(len(vx[:,i])):
                    peak = vx[k,i]
                    thisdx = dx[k,i]
                    dkmin = min(dkmin, thisdx)
                    if(thisdx == outdres): # this cell is unrefined
                        kincr = refinefac
                    elif(thisdx == outdres / 2):
                        kincr = int(refinefac / 2)
                    elif(thisdx == outdres / 4):
                        kincr = int(refinefac / 4)
                    else:
                        kincr = 1
                    # calculate the cumulative distribution of this line at each velocity bin edge
                    cdfs = 0.5 * (1 + special.erf((binvals - peak) / erfdenom)) * weight[k,i] * kincr
                    #erfvals[:] = [closest_erf_value(erfx, erfy, vval) for vval in (binvals - peak) / erfdenom]
                    #cdfs = 0.5 * (1 + erfvals * weight[k,i])
                    # subtract adjacent values to get the contribution to each bin
                    hist = hist + np.diff(cdfs)
                    k += kincr
                    if(k == len(vx[:,i])):
                        break
                if(dkmin == outdres):
                    iincr = refinefac
                elif(dkmin == outdres / 2):
                    iincr = int(refinefac / 2)
                elif(dkmin == outdres / 4):
                    iincr = int(refinefac / 4)
                else:
                    iincr = 1
                # this next bit handles binning together a refinefac**2 patch into one output cell
                # jj==0 handles glomming together the direction perpindicular to the slices
                # i//refinefac ==0 handles glomming to gether along the slice
                #print jj, i//refinefac, jj==0,i//refinefac==0
                if(jj == 0 and i%refinefac == 0):    
                    thesehists.append(hist * iincr)
                else:
                    thesehists[i//refinefac] += hist * iincr
               # print 'incrimenting i by ',iincr
                i += iincr
                if(i == inres):
                    break
            # figure out if we can skip reading some of these slices
            if(mindx == outdres): # there are no refined cells in this slice
                # all the subslices are going to be the same
                jincr = refinefac
            elif (mindx == outdres / 2): # there is only one level of refinement in this slice
                # the first refinefac/2 subsclices are going to be the same
                jincr = int(refinefac / 2)
            elif (mindx == outdres / 4): # there are two levels of refinement in this slice
                # the first refinefac/4 subsclices are going to be the same
                jincr = int(refinefac / 4)
            else:
                jincr = 1
            thesehists *= jincr
            if(jj == 0):
                thesehistsaccum = copy.deepcopy(thesehists)
            else:
                thesehistsaccum += thesehists
            #print 'incrimenting j by ',jincr
            jj += jincr
            if(jj == refinefac):
                break;
            
    # once we have the histograms of mass-weighted velocity along each point for this
    # row, save it to an hdf5 file
    f = h5py.File(specdir+'spectra_N2Hplus_'+str(j).zfill(4)+'.hdf5', 'w')
    dset = f.create_dataset('spectraN2Hplus', data = thesehistsaccum, dtype='float32')
    dset.attrs['slowindex'] = j
    dset.attrs[sliceax] = outpty
    f.close()
    
    del(slc)
    del(frb)
    del(f)
    del(dset)
    del(x)
    del(vx)
    del(dx)
    del(rhoN2Hplus)
    del(weight)
    del(hist)
    del(thesehists)
    del(thesehistsaccum)
    gc.collect()
    
    
