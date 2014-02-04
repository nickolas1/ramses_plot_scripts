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
from os.path import expanduser
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

tc = '0.5'
tc1 = '0.9'

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)

fig = plt.figure(figsize = (8, 7))
ax0 = fig.add_axes([0.1, 0.1, 0.75, 0.25])
ax1 = fig.add_axes([0.1, 0.3635, 0.75, 0.25])
ax2 = fig.add_axes([0.1, 0.615, 0.75, 0.25])

#snap = int(sys.argv[1])
#axis = int(sys.argv[2])
#filnumber = int(sys.argv[3])
# hardwire this
snap = 81
axis = 0
filnumber = 0

fildir = 'reduced_'+str(snap).zfill(5)+'/filaments'+str(axis)+'_'+str(filnumber)+'/'

prefix = 'reduced_'+str(snap).zfill(5)+'/'
infoname = prefix+'info_'+str(snap).zfill(5)+'.txt'
filprefix = prefix+'filament_'+str(filnumber)+'/'
filname = prefix+'filaments'+str(axis)+'_'+str(snap).zfill(5)+'.txt'

(boxlen, unit_l) = get_boxsize(infoname)
if boxlen > 7:
    sdoff = np.log10(4)
    vdoff = np.log10(8)
else:
    sdoff = 0.0
    vdoff = 0.0
    
imshowmap = 'nickmapSD'
#imshowmap = 'bone_r'
cdmin = -4.11 - sdoff + 1.2
cdmax = 0.069 - sdoff 

(lmin, lmax) = get_level_min_max(infoname)
(boxlen, unit_l) = get_boxsize(infoname)

plotsurfacedensity = True

# read in the rectangles that define the filaments we're interested in
# these are in units of pixels in the finder image, so we will need to translate these
# to unitary units!
rectdata = ascii.read(filname)
fil = rectdata[filnumber]

leftpoint = np.array([fil[1], fil[2]])
rightpoint = np.array([fil[3], fil[4]])
width = fil[5]
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
print length

# we will move along the lower-left line of the rectangle, starting from startpoint 
# and ending at endpoint 
startpoint = leftpoint - orthovec * width/2
endpoint = rightpoint - orthovec * width/2
print startpoint,endpoint
# for the distance to move in each step, choose a step size closest to the number of
# grid cells (at the finest refinement level) it would take to traverse the box 
expandfac = 4 # resample the surface density map
nl = int(length * expandfac)
nw = int(width * expandfac)
dl = length / nl
dw = width / nw

print nl, nw

if plotsurfacedensity:
    subbox = np.zeros([nw, nl])

    file = prefix+'surface_density_'+str(axis)+'.hdf5'
    print snap,file
    f = h5py.File(file, 'r')
    sd = f['surface_density']

    ny = sd.shape[1]
    for il in xrange(nl):
        l = startpoint + vec * dl * (il + 0.5)
        for iw in xrange(nw):
            pt = l + orthovec * dw * (iw + 0.5)
            # the sd array and imshow have different row-column ordering 
            # convention than everything else in the world, including the 
            # coordinates that we use to define the rectangles. so the 
            # ordering of the points here is reversed.
            subbox[iw, il] = sd[int(pt[1]),int(pt[0])]
    f.close()
    to_unitary = 1 / 2**lmax
    lscale = length * to_unitary * 10 * boxlen  #unit length is 10 pc
    ax0.imshow(subbox,
        origin='lower',
        extent = [0.001, lscale, 0.001, lscale*width/length],
        vmin = cdmin,
        vmax = cdmax,
        cmap = imshowmap,
        interpolation = 'nearest')
    
    # turn off axes
    #ax0.set_frame_on(False)
    #ax0.axes.get_yaxis().set_visible(False)
    #ax0.axes.get_xaxis().set_visible(False)
    ax0.get_yaxis().tick_left()
    for line in ax0.xaxis.get_ticklines():
        line.set_color(tc1)
    for line in ax0.yaxis.get_ticklines():
        line.set_color(tc1) 
    for line in ax0.yaxis.get_ticklines():
        line.set_color(tc1)     
    ax0.tick_params(which = 'minor', color = tc1) 
    ax0.grid(False)
    
    ax0.set_xlabel(r'$\mathdefault{L_{fil}}$ / pc', fontproperties = tfm, size = 15, color=tc)
    ax0.set_ylabel(r'$\mathdefault{W_{fil}}$ / pc', fontproperties = tfm, size = 15, color=tc)
    for label in ax0.get_xticklabels() + ax0.get_yticklabels():
        label.set_fontproperties(tfm)
        label.set_color(tc)

"""
this section plots the 'observed' density-weighted line of sight velocity.
this is plotted in the middle axis.
"""
spectra = []
#for i in xrange(358):
for filename in glob.glob(fildir+'combinedspectrum_*.hdf5'):    
#    filename = 'combinedspectrum_'+str(i).zfill(4)+'.hdf5'
    print filename
    f = h5py.File(filename, 'r')
    spectrum = np.array(f['spectra'])
    # normalize each slice by the rms of itself. this allows you to cut by 'detection' 
    # rather than the 'intensity' of the line. comment this out to go back to intensity.
    # note: this is probably bogus
    #spectrum /= np.sqrt(np.mean(spectrum**2))
    spectra.append(spectrum)
    f.close()

spectra /= np.mean(spectra)
#spectra /= np.sqrt(np.mean(spectrum**2))

print np.min(spectra),np.mean(spectra),np.median(spectra),np.max(spectra)

f = h5py.File(fildir+'spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

ax1.imshow(np.transpose(spectra),
    interpolation='nearest',
    origin = 'lower',
    extent = [0.001, lscale, 0.001, np.max(vels) - np.min(vels)],
    aspect = 0.13,
    #vmin = 3,
    vmax = 8,
    #cmap = 'gray_r',
    cmap = 'nickmapVD')

ax1.get_yaxis().tick_left()
for line in ax1.xaxis.get_ticklines():
    line.set_color(tc1)
for line in ax1.yaxis.get_ticklines():
    line.set_color(tc1) 
for line in ax1.yaxis.get_ticklines():
    line.set_color(tc1)     
ax1.tick_params(which = 'minor', color = tc1) 
ax1.tick_params(labelbottom='off')
ax1.grid(color='0.0',alpha=0.1)


ax1.set_ylabel(r'v / km $\mathdefault{s^{-1}}$', fontproperties = tfm, size = 15, color=tc)
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontproperties(tfm)
    label.set_color(tc)    
    
""""""""""""""""""""""""""""
"""


f = h5py.File(fildir+'spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

l_det = [] # detected l coordinates (length)
w_det = [] # detected w coordinates (width)
v_det = [] # detected v coordinates (velocity)
#for i in xrange(358):

l = 0 # keep track of l coordinate
for filename in glob.glob(fildir+'spectra_*.hdf5'):    
#    filename = 'combinedspectrum_'+str(i).zfill(4)+'.hdf5'
    print filename
    f = h5py.File(filename, 'r')
    specs = f['spectra']
    # treat the spectrum of each point individually
    spectrum = np.zeros(specs.shape[1])
    
    for s in xrange(specs.shape[0]):
        spec = np.array(specs[s])
        # normalize the spectrum by its rms
        spec /= np.sqrt(np.mean(spec**2))
        
        # set 'detected' vels equal to one, the rest to zero
        spec[spec < 3] = 0.0
        spec[spec > 0] = 1.0
        
        # march through and glom together touching detections
        runningn = 0
        runningv = 0
        if spec[0] >= 1:
            runningn = spec[0]
            runningv = vels[0] * spec[0]
        for i in xrange(1,len(spec)):
            if (spec[i] >= 1) & (spec[i-1] == 0):
                # we are starting a new detection
                runningn = spec[i]
                runningv = vels[i] * spec[i]
            if (spec[i] >= 1) & (spec[i-1] >= 1):
                # we are continuing a new detection
                runningn += spec[i]
                runningv += vels[i] * spec[i]
            if (spec[i] == 0) & (spec[i-1] >= 1):
                # we are ending a detection, record it
                l_det.append(l*dl)
                w_det.append(s*dw)
                # need to change w to be an average.
                v_det.append(runningv / runningn)
                # check to see if this is really averaging- if runningn>1, print
                #if runningn > 1:
                #    print runningv, runningv/runningn, runningn
    l += 1
    f.close()
    
# v_det is already in km/s. 
# offset by vmin to match the middle plot
v_det -= np.min(vels)

#convert l_det and w_det to pc
to_unitary = 1 / 2**lmax * expandfac
lscale2 = length * to_unitary * 10 * boxlen / l  #unit length is 10 pc
l_det = np.array(l_det)*lscale2
w_det = np.array(w_det)*lscale2
#ax2.scatter(l_det,v_det,marker='.',s=13,facecolor='0.0',lw=0,alpha=.33)
#ax2.set_xlim(0.001, lscale)
#ax2.set_ylim(np.min(v_det)-.1, np.max(v_det)+.1)
#ax2.set_axis_bgcolor('1.0')

dvel = vels[1] - vels[0]
velmin = 2.1
velmax = 5.6
velbins = np.arange(velmin, velmax, dvel)
dl = (np.max(l_det) - np.min(l_det))/l
lbins = np.arange(np.min(l_det), np.max(l_det), dl)

ax2.tick_params(labelbottom='off')


H, xedges, yedges = np.histogram2d(v_det, l_det, bins=(velbins, lbins))
print np.min(H), np.max(H), np.mean(H), np.median(H)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
ax2.imshow(H,
    extent = [0.001, lscale, velmin, velmax],
    origin = 'lower',
    interpolation = 'nearest',
    aspect=.25,
    cmap = 'nickmapVD',
    vmax = 8
    )

ax2.get_yaxis().tick_left()
ax2.get_xaxis().tick_bottom()
for line in ax2.xaxis.get_ticklines():
    line.set_color(tc1)
for line in ax2.yaxis.get_ticklines():
    line.set_color(tc1) 
for line in ax2.yaxis.get_ticklines():
    line.set_color(tc1)     
ax2.tick_params(which = 'minor', color = tc1) 
ax2.tick_params(labelbottom='off')
ax2.grid(color='0.0',alpha=0.1)

ax2.set_ylabel(r'v / km $\mathdefault{s^{-1}}$', fontproperties = tfm, size = 15, color=tc)
for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_fontproperties(tfm)
    label.set_color(tc)

plt.savefig('foofil1.pdf')       
sys.exit()

plt.savefig('fooppv.png',dpi=200)
plt.savefig('fooppv.pdf')
#plt.show()