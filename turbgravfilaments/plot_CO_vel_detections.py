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

fig = plt.figure(figsize = (5, 3))
ax = fig.add_axes([.2, .2, .75, .75])

snap = int(sys.argv[1])
axis = int(sys.argv[2])
fileprefix = 'reduced_'+str(snap).zfill(5)+'/'


# read in the rectangle from the filament definition
rectdata = ascii.read(fileprefix+'filaments'+str(axis)+'_'+str(snap).zfill(5)+'.txt')
for fil in rectdata:
    leftpoint = np.array([fil[1], fil[2]])
    rightpoint = np.array([fil[3], fil[4]])
    width = fil[5]
    vector = rightpoint - leftpoint
    orthovec = (-vector[1], vector[0])
    orthovec /= np.linalg.norm(orthovec)
    x = (leftpoint[0], rightpoint[0])
    y = (leftpoint[1], rightpoint[1])
    ul = leftpoint + orthovec * width/2
    ll = leftpoint - orthovec * width/2
    ur = rightpoint + orthovec * width/2
    lr = rightpoint - orthovec * width/2
    rectangle = np.transpose([ul, ll, lr, ur, ul])

    # this rectangle is in unitary units.
    
    print ll
    print ur
    
    # read in the detections file
    # NOTE this is [z, y, [vels]]
    f = h5py.File(fileprefix+'posvel_'+str(axis)+'/detections.hdf5')
    dets = np.array(f['veldetections'])
    f.close()
    res = dets.shape[0]
    dres = 1 / res
    
    # get index limits of the box
    zlo = int(np.floor(ll[1] / dres))
    zhi = int(np.ceil(ur[1] / dres))
    ylo = int(np.floor(ll[0] / dres))
    yhi = int(np.ceil(ur[0] / dres))
   
    # march along y and collect detections over z
    ys = []
    collecteddets = []
    for y in xrange(ylo, yhi+1, 1):
        zsumdets = (dets[zlo:zhi+1, y]).flatten()
        zsumdets = zsumdets[np.logical_not(np.isnan(zsumdets))]
        collecteddets.append(zsumdets.tolist())
        ys.append((np.ones(len(zsumdets)) * (y + 0.5) * dres).tolist())
     
ys = np.array(reduce(lambda x,y: x+y,ys))
collecteddets = np.array(reduce(lambda x,y: x+y,collecteddets))

# plot a 2D histogram of the detected vels as a function of y
f = h5py.File(fileprefix+'posvel_'+str(axis)+'/spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

dvel = vels[1] - vels[0]
dvel *= 1
velmin = -2
velmax = 2
velbins = np.arange(velmin, velmax, dvel)
print len(ys)
print np.min(ys)
print (yhi+1 - ylo)
dy = (np.max(ys) - np.min(ys))/(yhi+1 - ylo)
ybins = np.arange(np.min(ys), np.max(ys), dy)

#ax.tick_params(labelbottom='off')

H, xedges, yedges = np.histogram2d(collecteddets, ys, bins=(velbins, ybins))
print np.min(H), np.max(H), np.mean(H), np.median(H)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
ax.imshow(H,
    extent = [leftpoint[0], rightpoint[1], velmin, velmax],
    origin = 'lower',
    interpolation = 'nearest',
    aspect = .03,
    cmap = 'gray_r',
    vmax = 30
    )


ax.get_yaxis().tick_left()
ax.get_xaxis().tick_bottom()
for line in ax.xaxis.get_ticklines():
    line.set_color(tc1)
for line in ax.yaxis.get_ticklines():
    line.set_color(tc1) 
for line in ax.yaxis.get_ticklines():
    line.set_color(tc1)     
ax.tick_params(which = 'minor', color = tc1) 
ax.tick_params(labelbottom='off')
ax.grid(color='0.0',alpha=0.1)

ax.set_ylabel(r'$\mathdefault{v_{los}}$ / km $\mathdefault{s^{-1}}$', fontproperties = tfm, size = 15, color=tc)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)
    label.set_color(tc)

plt.savefig('fil0.pdf')       
sys.exit()

    
    
    
    

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
#ax.scatter(l_det,v_det,marker='.',s=13,facecolor='0.0',lw=0,alpha=.33)
#ax.set_xlim(0.001, lscale)
#ax.set_ylim(np.min(v_det)-.1, np.max(v_det)+.1)
#ax.set_axis_bgcolor('1.0')

dvel = vels[1] - vels[0]
velmin = 2.1
velmax = 5.6
velbins = np.arange(velmin, velmax, dvel)
dl = (np.max(l_det) - np.min(l_det))/l
lbins = np.arange(np.min(l_det), np.max(l_det), dl)

ax.tick_params(labelbottom='off')


H, xedges, yedges = np.histogram2d(v_det, l_det, bins=(velbins, lbins))
print np.min(H), np.max(H), np.mean(H), np.median(H)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
ax.imshow(H,
    extent = [0.001, lscale, velmin, velmax],
    origin = 'lower',
    interpolation = 'nearest',
    aspect=.25,
    cmap = 'gray_r',
    vmax = 8
    )

ax.get_yaxis().tick_left()
ax.get_xaxis().tick_bottom()
for line in ax.xaxis.get_ticklines():
    line.set_color(tc1)
for line in ax.yaxis.get_ticklines():
    line.set_color(tc1) 
for line in ax.yaxis.get_ticklines():
    line.set_color(tc1)     
ax.tick_params(which = 'minor', color = tc1) 
ax.tick_params(labelbottom='off')
ax.grid(color='0.0',alpha=0.1)

ax.set_ylabel(r'$\mathdefault{v_{los}}$ / km $\mathdefault{s^{-1}}$', fontproperties = tfm, size = 15, color=tc)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)
    label.set_color(tc)

plt.savefig('fil0.pdf')       
sys.exit()

plt.savefig('fooppv.png',dpi=200)
plt.savefig('fooppv.pdf')
#plt.show()