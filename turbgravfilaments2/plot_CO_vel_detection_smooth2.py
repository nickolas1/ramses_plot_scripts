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
import glob
from astropy.io import ascii
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from scipy import special
from mpl_toolkits.mplot3d import Axes3D

"""
usage:
python fooppv.py N A F
N: number of output to use. reduced_N needs to be here.
A: axis of the projection (0, 1, 2)
F: filament number 
"""

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

mpl.rcParams['xtick.major.size'] = 9
mpl.rcParams['axes.unicode_minus'] = False

tc = '0.4'
tc1 = '0.9'

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)

fig = plt.figure(figsize = (7, 3))

snap = int(sys.argv[1])
axis = int(sys.argv[2])
fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

infoname = fileprefix+'info_'+str(snap).zfill(5)+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)
print boxlen, unit_l
unitarytopc = boxlen * unit_l / const.pc.cgs.value

# read in the rectangle from the filament definition
rectdata = ascii.read(fileprefix+'filaments'+str(axis)+'_'+str(snap).zfill(5)+'.txt')
for fil in rectdata:
    filnum = fil[0]
    leftpoint = np.array([fil[1], fil[2]])
    rightpoint = np.array([fil[3], fil[4]]) 
    width = fil[5] 
    vector = rightpoint - leftpoint
    length = np.linalg.norm(vector)
    orthovec = (-vector[1], vector[0])
    orthovec /= np.linalg.norm(orthovec)
    vector /= np.linalg.norm(vector)
    x = (leftpoint[0], rightpoint[0])
    y = (leftpoint[1], rightpoint[1])
    ul = leftpoint + orthovec * width/2
    ll = leftpoint - orthovec * width/2
    ur = rightpoint + orthovec * width/2
    lr = rightpoint - orthovec * width/2
    rectangle = np.transpose([ul, ll, lr, ur, ul])
    # this rectangle is in unitary units.
    
    velmin = 1.0 * fil[6]
    velmax = 1.0 * fil[7]
    veltickspace = 1.0 * fil[8]
    print velmin, velmax, veltickspace
    
    print 'lower left corner: ',ll
    print 'upper right corner: ',ur

    # we will move along the lower-left line of the rectangle, starting from startpoint 
    # and ending at endpoint 
    startpoint = leftpoint - orthovec * width/2
    endpoint = rightpoint - orthovec * width/2
    print 'startpoint, endpoint: ',startpoint,endpoint
    
    # the box size is 10 pc. convert the length to physical scale
    lscale = length * unitarytopc
    
    if length * unitarytopc < 3:
        xtickspan = 0.5
    else:
        xtickspan = 1.0
    xtickvals = np.arange(0, length * unitarytopc, xtickspan)

    """
    this section makes a resampled, rotated column density map of the filament
    that we're dealing with.
    """
    plotSurfaceDensity = False
    if plotSurfaceDensity:
        ax0 = fig.add_axes([0.15, 0.1, 0.8, 0.8])
        # read in the surface density file and figure out the spatial resolution
        file = fileprefix+'surface_density_CO'+str(axis)+'.hdf5'
        print snap,file
        f = h5py.File(file, 'r')
        sd = f['surface_density_CO']

        sd = 10**np.array(sd)  # convert to linear units
        print np.mean(sd),np.max(sd)
        sd /= (2.33 * const.m_p.cgs.value) # convert to number density
        print np.mean(sd),np.max(sd)
        sd /= (2.42e14 * 2.94e6) # non-temperature factors of IC18O conversion
        sd /= (np.exp(5.27/10) / (np.exp(5.27/10) - 1)) # temperature part
        print np.mean(sd),np.max(sd)
        cdmin = 0
        cdmax = 5

        res_sd = sd.shape[1]
        dres_sd = 1 / res_sd
    
        expandfac = 12 # resample the surface density map to make off-axis pixels look nice
        nl = int(length * res_sd * expandfac)
        nw = int(width * res_sd * expandfac)
        dl = length / nl
        dw = width / nw
        print 'nl, nw: ',nl, nw
    
        # set up empty array for the resampled surface density map
        subbox = np.zeros([nw, nl])
    
        for il in xrange(nl):
            l = startpoint + vector * dl * (il + 0.5)
            for iw in xrange(nw):
                pt = l + orthovec * dw * (iw + 0.5)
                pt /= dres_sd
                # the sd array and imshow have different row-column ordering 
                # convention than everything else in the world, including the 
                # coordinates that we use to define the rectangles. so the 
                # ordering of the points here is reversed.
                subbox[iw, il] = sd[int(np.remainder(pt[1], res_sd)),int(np.remainder(pt[0], res_sd))]
        f.close() # close the surface density file
        print np.min(subbox), np.max(subbox)
    
        imshowmap = 'nickmapVD2'
        ax0.imshow(subbox,
            origin='lower',
            extent = [0, lscale, 0, lscale*width/length],
            vmin = cdmin,
            vmax = cdmax,
            cmap = imshowmap,
            interpolation = 'nearest')
    
        if length * unitarytopc < 3:
            xtickspan = 0.5
        else:
            xtickspan = 1.0
        ax0.xaxis.set_ticks(xtickvals)
        ax0.yaxis.set_ticks(np.arange(0, width * unitarytopc, 0.2))

        ax0.get_yaxis().tick_left()
        for line in ax0.xaxis.get_ticklines():
            line.set_color(tc1)
        for line in ax0.yaxis.get_ticklines():
            line.set_color(tc1) 
        for line in ax0.yaxis.get_ticklines():
            line.set_color(tc1)     
        ax0.tick_params(which = 'minor', color = tc1) 
        ax0.grid(False)
    
        ax0.set_xlabel(r'L / pc', fontproperties = tfm, size = 15, color=tc)
        ax0.set_ylabel(r'W / pc', fontproperties = tfm, size = 15, color=tc)
        for label in ax0.get_xticklabels() + ax0.get_yticklabels():
            label.set_fontproperties(tfm)
            label.set_color(tc)
        
        plt.savefig('filsd'+str(filnum)+'.png',dpi=400) 
        plt.savefig('filsd'+str(filnum)+'.pdf',dpi=400) 
        plt.clf()
    
    """
    this section makes a plot of the density-weighted line of site velocity
    along the filament. the velocity is summed along the W filament direction.
    """
    ax0 = fig.add_axes([0.15, 0.04, 0.8, 0.8])
    f = h5py.File(fileprefix+'posvel_'+str(axis)+'/spectrumvels.hdf5')
    vels = np.array(f['binmidskms'])
    f.close()
    print 'length of spectra: ',len(vels)
    specfile = fileprefix+'posvel_'+str(axis)+'/spectra_0000.hdf5'
    f = h5py.File(specfile)
    specs = f['spectra']
    res_spec = specs.shape[0]
    dres_spec = 1 / res_spec
    f.close()
    
    # read in the detections file and figure out the spatial resolution of the spectra data
    # NOTE this is [z, y, [vels]]
    #f = h5py.File(fileprefix+'posvel_'+str(axis)+'/detections.hdf5')
    #dets = np.array(f['veldetections'])
    #f.close()
    #res_spec = dets.shape[0]
    #dres_spec = 1 / res_spec
    
    # get index limits of the box
    zlo = int(np.floor(min(ll[1], lr[1]) / dres_spec))
    zhi = int(np.ceil(max(ul[1], ur[1]) / dres_spec))
    ylo = int(np.floor(min(ll[0], ul[0]) / dres_spec))
    yhi = int(np.ceil(max(lr[0], ur[0]) / dres_spec))
    spany = yhi - ylo
    spanz = zhi - zlo
    print 'untransformed box spany, spanz: ',spany, spanz
    print zlo, zhi
    
    
    # for the distance to move in each step, choose a step size closest to the number of
    # grid cells (at the finest refinement level) it would take to traverse the box 
    expandfac = 1 # resample the surface density map
    nl = int(length * res_spec * expandfac)
    nw = int(width * res_spec * expandfac)
    dl = length / nl
    dw = width / nw  
    print 'nl, nw: ',nl, nw
 
 
    # set up the array for the untransformed spectra
    specfield = np.zeros([spany, spanz, len(vels)]) 
    
    # get all of the spectra in this box into an array for easy access
    for i in xrange(zlo, zhi, 1):
        specfile = fileprefix+'posvel_'+str(axis)+'/spectra_'+str(np.remainder(i, res_spec)).zfill(4)+'.hdf5'
        f = h5py.File(specfile)
        specs = f['spectra']
        
        for s in xrange(ylo, yhi, 1):
            spec = np.array(specs[s])
            specsmooth = np.zeros(len(vels))
            dvel = vels[1] - vels[0]
            veledges = np.arange(vels[0] - dvel/2, vels[-1] + dvel, dvel)
            print zlo, i, zhi,"    ",ylo, s, yhi
            for v in xrange(len(vels)):
                peak = vels[v]
                sigma = 0.052
                amp =spec[v] 
                if amp > 1.e-100:
                    erfdenom = np.sqrt(2*sigma**2)
                    cdfs = 0.5 * (1 + special.erf((veledges - peak) / erfdenom)) * amp
                    weights = np.roll(cdfs, -1) - cdfs
                    specsmooth = specsmooth + weights[:-1]         
            specfield[s-ylo, i-zlo] = specsmooth   
        f.close() 
    
    # set up an empty array of spectra along the length of the filament
    transspec = np.zeros([nl, len(vels)])
    
    print 'hello'
    # march along the length and collect spectra along width
    for il in xrange(nl):
        l = startpoint + vector * dl * (il + 0.5)
        for iw in xrange(nw):
            pt = l + orthovec * dw * (iw + 0.5)
            # pt is now a y-z position. convert this to coordinates 
            pt /= dres_spec
            iy = int(pt[0])
            iz = int(pt[1])
            #print specfield[iy-ylo][iz-zlo]
            transspec[il,:] += specfield[iy-ylo][iz-zlo]
          #  print specfield[iy-ylo][iz-zlo]
        #print transspec[il,:]
    # set up the length and velocity bins, and make the histogram
    # use length / dres to make length bins. i.e. the subsampling 
    # used to march along l only serves to weight each pixel appropriately
    # in the transformed space- plot at the same resolution as the simulation.
    
    dvel = vels[1] - vels[0]
    dvel *= 1
   # velmin = -2
   # velmax = 2
    transspecsub = transspec[:,(vels >= velmin) & (vels <= velmax)] 


    transspecsub /= np.max(transspecsub)
    print np.min(transspecsub), np.max(transspecsub), np.median(transspecsub), np.mean(transspecsub)
    transspecsub_nonzero = transspecsub[transspecsub > 0.01]
    print np.min(transspecsub_nonzero), np.median(transspecsub_nonzero), np.max(transspecsub_nonzero)
    if(snap == 29): # reverse L dimension for this filament
        if(filnum == 0):
            transspecsub = np.swapaxes(np.swapaxes(transspecsub,0,0)[::-1],0,0)
    ax0.imshow(transspecsub.T,
        extent = [0, lscale, velmin, velmax],
        origin = 'lower',
        interpolation = 'nearest',
        aspect = 'auto',
        cmap = 'nickmapSD2',
        vmin = 0.01,
        vmax = min(1.0, 5. * np.median(transspecsub_nonzero))
       # cmap = 'gray_r',
       # vmin = -1,
       # vmax = 0
        )
    ax0.xaxis.set_ticks(xtickvals)
    velticksvals = np.union1d(np.arange(0, velmin, -veltickspace), np.arange(0, velmax, veltickspace))
    ax0.yaxis.set_ticks(velticksvals)
    
    ax0.get_yaxis().tick_left()
    ax0.get_xaxis().tick_bottom()
    for line in ax0.xaxis.get_ticklines():
        line.set_color(tc1)
    for line in ax0.yaxis.get_ticklines():
        line.set_color(tc1) 
    for line in ax0.yaxis.get_ticklines():
        line.set_color(tc1)     
    ax0.tick_params(which = 'minor', color = tc1) 
    ax0.tick_params(labelbottom='off')
    ax0.grid(color='0.0',alpha=0.1)
    
    ax0.set_ylabel(r'$\mathdefault{v_{los}}$ / km $\mathdefault{s^{-1}}$', fontproperties = tfm, size = 15, color=tc)
    for label in ax0.get_xticklabels() + ax0.get_yticklabels():
        label.set_fontproperties(tfm)
        label.set_color(tc)
    
    #plt.savefig('filvel'+str(filnum)+'.pdf')
    
    # now make a rasterized version with no plot labels
    for label in ax0.get_xticklabels() + ax0.get_yticklabels():
        label.set_color('white')
    ax0.set_ylabel(r'', fontproperties = tfm, size = 15, color='white', zorder=50)
    ax0.set_rasterized(True)
    plt.savefig('filvelsmooth'+str(filnum)+'.png', dpi=400)
    plt.clf()
    

   