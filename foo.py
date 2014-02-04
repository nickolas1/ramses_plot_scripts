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

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)


fig = plt.figure(figsize = (5, 3.5))
ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])

snap = int(sys.argv[1])
axis = int(sys.argv[2])
filnumber = int(sys.argv[3])-1

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
cdmin = -4.11 - sdoff
cdmax = 0.069 - sdoff

(lmin, lmax) = get_level_min_max(infoname)
(boxlen, unit_l) = get_boxsize(infoname)

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

ax.imshow(subbox,
    origin='lower',
    vmin = cdmin,
    vmax = cdmax,
    cmap = imshowmap)
    
# turn off axes
ax.set_frame_on(False)
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)

plt.show()
        
sys.exit()

spectra = []
#for i in xrange(358):
for filename in glob.glob('combinedspectrum_*.hdf5'):    
#    filename = 'combinedspectrum_'+str(i).zfill(4)+'.hdf5'
    f = h5py.File(filename, 'r')
    spectrum = np.array(f['spectrum'])
    spectra.append(spectrum)
    f.close()

spectra /= np.median(spectra)
    
f = h5py.File('spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

ax.imshow(np.transpose(spectra),
    interpolation='nearest',
    origin = 'lower',
    extent = [0, 1, np.min(vels), np.max(vels)],
    aspect = 0.1,
    vmin = 3,
    vmax = 6,
    cmap = 'gray_r')

print np.min(spectra)
print np.max(spectra)
print np.mean(spectra)
print np.median(spectra)    

plt.savefig('foo.png')
#plt.show()