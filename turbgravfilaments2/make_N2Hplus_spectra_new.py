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


"""
    read in the surface density map. we will use this to normalize the spectra
    so that they are in approximately physical units
    
    convert surface density to integrated line intensity using
    Alvaro's empirical conversion:
    
    N(H2) [cm-2] ~ 7E21 * I(N2H+) [K kms-1]
"""
file = 'reduced_'+str(snap).zfill(5)+'/surface_density_C18O2.hdf5'
f = h5py.File(file, 'r')
sd = f['surface_density_C18O']
sd = 10**np.array(sd)  # convert to linear units
sd /= (2.33 * 1.672621777e-24) # convert to number density
sd /= 7.0e21
f.close()


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



vmax = 2.5e5 + 7.0e5
vmin = -2.5e5 - 8.0e5
# roughly match hacar et al by takin 0.05 km/s bins
bins = (vmax - vmin) / 1.e5 / 0.05
binvals = np.arange(vmin, 1.000001*vmax, (vmax - vmin) / bins)
binmids = 0.5 * (np.roll(binvals, -1) + binvals)
binmids = binmids[:len(binmids) - 1]
# get a version of the bins in km/s instead of cgs
binmidskms = binmids / 1.e5
    
# save the velocities to a file
f = h5py.File(specdir+'spectrumvelsN2Hplus.hdf5', 'w')
dset = f.create_dataset('binmidskms', data = binmidskms, dtype='float32')
f.close()

"""
    to keep things manageable, make this output map on the 1024**3 base grid.
    outres: the output resolution.
    inres: the resolution we sample the grid at.
"""
outres = 2**lmin
outdres = 1.0 / outres

refinefac = 2**(lmax - lmin)
inres = outres * refinefac
indres = 1.0 / inres

sigmaC18O = 0.0535 # thermal width of N2Hplus line in km/s
sigma = sigmaC18O * 1.e5 # convert to cm/s
erfdenom = np.sqrt(2*sigma**2)

"""
    centroid velocities and relative intensities of the N2H+ hyperfine lines
""" 
N2HplusCenters = np.array([6.9360, 5.9842, 5.5452, 0.9560, 0.0, -0.6109, -8.0064])
N2HplusCenters *= 1.e5 # convert to cm/s
N2HplusIntensities = np.array([0.03704, 0.18519, 0.11111, 0.18519, 0.25926, 0.11111, 0.11111])

for sj in xrange(200):
    outpty = (sj + 0.5) * outdres
    thesehistsaccum = np.zeros([outres, len(binmids)])
    print sj, outpty
    
    if(sd[sj,:].sum() == 0):
        print 'skipping ',sj
        continue
    
    j = 0
    for ij in xrange(refinefac):
        thesehists = np.zeros([outres, len(binmids)])
        inpty = (sj*refinefac + j + 0.5) * indres
        print 'inpty: ',inpty
        # get a slice
        slc = ds.h.slice(sliceax, inpty)
    
        # get it into a frb
        frb = slc.to_frb(
            (1.0, 'unitary'),           # get the whole extent of the box
            inres,                       # don't degrade anything
            center = [0.5, 0.5, 0.5],   # centered in the box
            height = (1.0, 'unitary'))  # get the whole extent of the box
    
        weight = np.array(frb['N2Hplus'])
    
        x = np.array(frb[los])
        vx = np.array(frb[vlos])
        dx = np.array(frb[dlos])
        mindx = np.min(dx)
        print 'max(dx), min(dx), outdres: ',np.max(dx),np.min(dx),outdres
        print 'max(rho), min(rho), outdres: ',np.max(weight),np.min(weight),outdres

        i = 0
        for ii in xrange(inres):
            # for each point along the slice, march along the projecting dimension
            # and turn each detection into a gaussian. bin this gaussian into the 
            # velbins.
            hist = np.zeros(len(binmids))
            k = 0
            dkmin = np.min(dx[:,i])
            
            if(weight[:,i].sum() > 0):
             #   print 'whaoooooooo ',i, i//refinefac,weight[:,i].sum()
                for ik in xrange(len(vx[:,i])):
                    kincr = 1
                    if weight[ik,i] > 0:
                        peak = vx[ik,i]
                        thisdx = dx[ik,i]
                        kincr = int(thisdx / indres)
                        for iline in xrange(7):
                            hyppeak = peak + N2HplusCenters[iline]
                            hypintens = weight[ik,i] * N2HplusIntensities[iline]
                            # calculate the cumulative distribution of this line at each velocity bin edge
                            cdfs = 0.5 * (1 + special.erf((binvals - hyppeak) / erfdenom)) * hypintens * kincr
                            # subtract adjacent values to get the contribution to each bin
                            hist = hist + np.diff(cdfs)
                      #  print 'hist ',hist
                    k += kincr
                    if(k == len(vx[:,i])):
                        break
                        
            iincr = int(dkmin / indres)
            # this next bit handles binning together a refinefac**2 patch into one output cell
            # ii//refinefac == 0 handles glomming to gether along the slice
            thesehists[i//refinefac] += hist * iincr
           # print 'incrimenting i by ',iincr
           # print i, i//refinefac,dkmin / indres,inres
            i += iincr
            if(i == inres):
                break
                
        jincr = int(mindx / indres)
        thesehistsaccum += thesehists * jincr
        #print 'incrementing j by ',jincr
        j += jincr
        if(j == refinefac):
            break

    # normalize to put into K
    integratedI = thesehistsaccum.sum(axis = 1)
    for m in xrange(outres):
        if integratedI[m] == 0 and sd[sj, m] > 0:
            print 'trouble! ',m, foo[m], sd[sj, m]
    normalisations = sd[sj, :] / integratedI
    thesehistsaccum *= normalisations[:, np.newaxis]
    thesehistsaccum = np.nan_to_num(thesehistsaccum)
                
    # once we have the histograms of mass-weighted velocity along each point for this
    # row, save it to an hdf5 file
    f = h5py.File(specdir+'spectra_N2Hplus_'+str(sj).zfill(4)+'.hdf5', 'w')
    dset = f.create_dataset('spectraN2Hplus', data = thesehistsaccum, dtype='float32')
    dset.attrs['slowindex'] = sj
    dset.attrs[sliceax] = outpty
    f.close()
    
    del(slc)
    del(frb)
    del(f)
    del(dset)
    del(x)
    del(vx)
    del(dx)
    del(weight)
    del(hist)
    del(thesehists)
    del(thesehistsaccum)
    gc.collect()
