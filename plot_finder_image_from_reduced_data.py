from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.colors as col
import matplotlib.cm as cm
import gc
import sys
import h5py
from astropy.io import ascii
from os.path import expanduser


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# the limits on the surface density color map is set up for the compact clouds. here we
# see if we're looking at diffuse clouds, and if so we adjust for that.
snapstr = str(int(sys.argv[1])).zfill(5)
infoname = 'reduced_'+snapstr+'/info_'+snapstr+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)
if boxlen > 7:
    sdoff = np.log10(4)
    vdoff = np.log10(8)
else:
    sdoff = 0.0
    vdoff = 0.0

snap = int(sys.argv[1])
axis = int(sys.argv[2])

   
imshowmap = 'nickmapSD'
#imshowmap = 'bone_r'
cdmin = -4.11 - sdoff
cdmax = 0.069 - sdoff

fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.out'


file = fileprefix+'surface_density_'+str(axis)+'.hdf5'

print snap,file
f = h5py.File(file, 'r')
sd = f['surface_density']
fig = plt.figure(figsize = (sd.shape[0]/200, sd.shape[1]/200), dpi=200)
ax = fig.add_axes([0., 0., 1., 1.])
ax.imshow(sd,
      interpolation = 'nearest',
      origin = 'lower',
      vmin = cdmin,
      vmax = cdmax,
      cmap = imshowmap)
# turn off axes
ax.set_frame_on(False)
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)

# see if we have any sink particles to plot
try:
    with open(sinkname): 
        sinks = get_sinks(sinkname)
        if len(sinks) > 0:
            # figure out the size of the sinks in units of 0-1
            #mincell = 1.0/2**lmax
            #sinkrad = 1.5 * mincell
    
            sinkpos = sinks[:,2:5]
            sinkpos[:] /= boxlen # shrink to 0-1 in all dimensions
            # get projected positions
            keep = np.array([1,1,1])
            keep[axis] = 0
            keep = np.array(keep, dtype=bool)
            sinkpos = sinkpos[:,keep]

            # restrict to same region as density plot
            #ledge = cntr[0] - wd/2
            #bedge = cntr[2] - ht/2
            #sinkpos[:] -= np.array([ledge, bedge])
            # convert to imshow scale
            #sinkpos *= res[1] / wd
            #sinkrad *= res[1] / wd
            sinkpos *= sd.shape[0]
            print sinkpos      
            sinkmass = sinks[:,1]
            # color by the log of mass. the minimum that we plot is 0.1 Msun,
            # max is a few hundred.
            mmin = np.log10(1)
            mmax = np.log10(100)
            sinkmap = cm.get_cmap('nickmapSink')
            sinkcolors = sinkmap((np.log10(sinkmass) - mmin) / (mmax - mmin))     
            ax.autoscale(False)
            #for s in xrange(len(sinks)):
            #    ax.add_artist(Circle((sinkpos[s,0],sinkpos[s,1]),sinkrad,fc=csink))
            ax.scatter(sinkpos[:,0],sinkpos[:,1],marker='.',s=9,facecolor=sinkcolors,lw=0)              
except IOError:
    pass     

x = [1,10, 100,1000]
y = [10,10,10,10]
#ax.scatter(x, y, s=.5)  
# define a rectangle by drawing a line on a filament and choosing a width
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
    #ax.plot(x,y,lw=.3,color=cred)
    ax.plot(rectangle[0], rectangle[1], lw=.7,color='m', solid_joinstyle='miter')
    

framesdir = 'finderimage'+str(axis)+'/'
if not os.path.exists(framesdir):
    os.makedirs(framesdir)

framename = 'finderimage'+str(axis)+'_frame_'+str(snap).zfill(5)+'.pdf'
print framename
plt.savefig(framename, dpi = 200)
f.close() 
plt.close() 
del(f)
del(sd)
gc.collect()
        
