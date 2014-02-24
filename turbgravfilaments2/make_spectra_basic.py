from __future__ import division

from yt.mods import *
import numpy as np
import gc
import sys
import h5py
import shutil
from os.path import expanduser
from scipy import special

"""
usage: python make_spectra_basic.py <snap> <axis>
e.g. python make_spectra_basic.py 18 2 will get create spectra along the z axis for
output_00018

this requires h5py. if you don't want to mess around with hdf5, save the data
in some other format

important comments about how things work, or things you might want to change, are in this 
triple-quotes style
"""

"""
    these definitions are for looking at certain density ranges
"""
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
    
# grab the command line arguments
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

""" 
    create the directory you want to put the spectra files in
"""
specdir = 'reduced_'+str(snap).zfill(5)+'/posvel_'+str(axis)+'/'
if not os.path.exists(specdir):
    os.makedirs(specdir)

ds = load(infoname)

""" 
    comment both of these out if you aren't using a density range
"""
# add new density fields
add_field('C18O', function=_C18O)
#add_field('N2Hplus', function=_N2Hplus)

"""
    set the velocity limits of your spectra in cgs
"""
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
    refinefac: how far down the refinement hierarchy to go. this choice creates:
            inres: the resolution we sample the grid at.
"""
outres = 1024
refinefac = 8

outdres = 1.0 / outres
inres = outres * refinefac
indres = 1.0 / inres

"""
    start at the bottom of the projection (y = 0), and create slices through the box in 
    the middle of each grid cell. we'll create an output file for each slice. if we're 
    interested in z velocities, the slices will be at given values of y, and will be of 
    the x-z plane. 
    
    march along the x direction of the slice, and create a spectra for that (x,y) point
    by taking each value of v and rho along the z direction of the slice.
"""

for j in xrange(outres):
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
            inres,                      # don't degrade anything
            center = [0.5, 0.5, 0.5],   # centered in the box
            height = (1.0, 'unitary'))  # get the whole extent of the box
    
        """
            if you just want density weighted velocity rather than a density subrange,
            replace this with np.array(frb['rho'])
            
            replace sigma with some sort of line width to smooth the velocities
        """
        rhoC18O = np.array(frb['C18O'])
       # rhoN2Hplus = np.array(frb['N2Hplus'])
        sigmaC18O = 0.0526 # thermal width of C18O line in km/s
    
        sigma = sigmaC18O * 1.e5 # convert to cm/s
        erfdenom = np.sqrt(2*sigma**2)  #denominator for the error function involved in 
                                        #calculating the cumulative distribution of a gaussian
        x = np.array(frb[los])      # the position along the line of sight
        vx = np.array(frb[vlos])    # the line of sight velocities
        dx = np.array(frb[dlos])    # the size of the cell at each point
        mindx = np.min(dx)          # the minimum cell size in this slice
        print 'max(dx), min(dx), outdres: ',np.max(dx),np.min(dx),outdres
        
        weight = indres * rhoC18O  
        # we need to grab rows from the slice differently depending on what axis we're projecting
        hist = np.zeros(len(binmids))
        erfvals = np.zeros(len(binvals))
        """
            this if axis === 0 part is very wrong .... 
        """
        if axis == 0:  
            for i in xrange(inres):
                hist, binedges = np.histogram(
                    vx[i,:],
                    range = (vmin, vmax),
                    bins = binvals,
                    weights = weight[i,:])
                thesehists.append(hist)
        """
            this if axis > 0 part is up to date
        """
        if axis > 0:
            i = 0
            for ii in xrange(inres):
                # for each point along the slice, march along the projecting dimension
                # and turn each detection into a gaussian. bin this gaussian into the 
                # velbins.
                hist[:] = 0
                k = 0
                dkmin = 1 # keep track of the minimum cell size along this line of sight
                if(weight[:,i].sum() >= 0): # if there isn't any interesting density in this slice, leave the spectra at zeros
                    for ik in xrange(len(vx[:,i])):
                        peak = vx[k,i]
                        thisdx = dx[k,i]
                        dkmin = min(dkmin, thisdx)
                        """
                            this if elif elif else thing lets you skip needless computations
                            if the cell we're looking at is unrefined. 
                            we do this in all three dimensions. this saves a great deal 
                            of time if the simulation is large
                        """
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
                """
                    this next bit handles binning together a refinefac**2 patch into one output cell
                    jj==0 handles glomming together the direction perpindicular to the slices
                    i%refinefac==0 handles glomming to gether along the slice
                """
                if(jj == 0 and i%refinefac == 0):    
                    thesehists.append(hist * iincr)
                else:
                    thesehists[i//refinefac] += hist * iincr
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
            if(jj == 0):
                thesehistsaccum = np.array(thesehists) * jincr
            else:
                thesehistsaccum += np.array(thesehists) * jincr
            jj += jincr
            if(jj == refinefac):
                break;
    """
        once we have the histograms of mass-weighted velocity along each point for this
        row, save it. I use hdf5 files, you can save it however you want
    """
    f = h5py.File(specdir+'spectra_C18O_'+str(j).zfill(4)+'.hdf5', 'w')
    dset = f.create_dataset('spectraC18O', data = thesehistsaccum, dtype='float32')
    dset.attrs['slowindex'] = j
    dset.attrs[sliceax] = outpty
    f.close()
     
    """
        force deletions for memory preservation
    """
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
    del(thesehistsaccum)
    gc.collect()
    
    
