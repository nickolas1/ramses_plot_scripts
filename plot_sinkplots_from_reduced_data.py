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
from os.path import expanduser
from matplotlib.collections import LineCollection

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)

# the limits on the surface density color map is set up for the compact clouds. here we
# see if we're looking at diffuse clouds, and if so we adjust for that.
infoname = 'reduced_00001/info_00001.txt'
(boxlen, unit_l) = get_boxsize(infoname)
if boxlen > 7:
    tlimmin = 3.5
    mlimmax = 1.e4
else:
    tlimmin = 0.5
    mlimmax = 1.e4


# these two from http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #segments2 = np.concatenate([0.999*points[:], 1.001*points[:]], axis=1)
    #segments = np.concatenate([segments1, segments2])
    
   # segments[:,0,0] *= 0.999
   # segments[:,1,0] *= 1.001
    
    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=6, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth,alpha=alpha)
    
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
 
          
"""
This second part makes plots of the sink particle masses with time
"""

# override the defaults for this movie plot
mpl.rc('grid', color='0.15')
mpl.rc('grid', linewidth='1.0')
mpl.rc('axes', facecolor='0.0')
mpl.rc('xtick', color='0.6')
mpl.rc('ytick', color='0.6')
mpl.rc('figure', facecolor='0.0')
mpl.rc('savefig', facecolor='0.0')

times = [0.0]
alltimes = []
sinkmasses = [0.0]

indivnames = []
indivmasses = []
indivtimes = []

snapstart = int(sys.argv[1])
snapend = int(sys.argv[2])
snapiter = int(sys.argv[3])

for snap in range(snapstart, snapend, snapiter):
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    framesdir = outdir+'sinkmasses/'
    if not os.path.exists(framesdir):
            os.makedirs(framesdir)
    
    infoname = fileprefix+'info_'+str(snap).zfill(5)+'.txt'
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.out'

    # see if we have any sink particles to plot
    try:
        with open(sinkname): 
            sinks = get_sinks(sinkname)
            (time, unit_t) = get_time(infoname)
            timeMyr = time * unit_t / 31557600.0 / 1.e6
            alltimes.append(timeMyr)
            if len(sinks) > 0:
                sinkmass = sinks[:,1].sum() # total mass of sinks in Msun
                times.append(timeMyr)
                sinkmasses.append(sinkmass)
                timerow = np.ones((sinks.shape[0],1)) * timeMyr
                # nmt = name, mass, time
                nmt = np.hstack((sinks[:,[0,1]],timerow))
                indivnames.append(nmt[:,0])
                indivmasses.append(nmt[:,1])
                indivtimes.append(nmt[:,2])
    except IOError:
        pass      

# flatten lists of individual sink properties, stick in an array
indivnames = [j for i in indivnames for j in i]
indivmasses = [j for i in indivmasses for j in i] 
indivtimes = [j for i in indivtimes for j in i]     
# nmt is name mass time
nmt = np.array([indivnames, indivmasses, indivtimes])
    
times = np.array(times[1:])
sinkmasses = np.array(sinkmasses[1:])    

print 'times = ',times
print 'total sink masses = ',sinkmasses

# plot individual sink masses
mmin = np.log10(1)
mmax = np.log10(100)

for snap in range(snapstart, snapend, snapiter):
    fig = plt.figure(figsize = (5, 3.5))
    ax = fig.add_axes([0.248, 0.2, 0.652, 0.75])
    t = alltimes[snap-1]
    print 'plotting sink masses for time ',t
    
    nplotted = 0
    for i in xrange(int(max(nmt[0,:]))):
        thisone = nmt[:,nmt[0,:] == i]
        if len(thisone[0]) == 0: continue
        x = thisone[2,:]
        if np.min(x) > t: continue
        nplotted += 1
        sel = (x <= t)
        x = x[sel]
        y = thisone[1,sel]
   # if len(x) < 50 and len(x) > 1:
   #     xnew = np.linspace(np.min(x), np.max(x), 50)
   #     ynew = np.interp(xnew, x, y)       
   #     x = xnew
   #     y = ynew
        colorline(x, y, 
            np.log10(y), 
            cmap='nickmapSink', 
            norm=plt.Normalize(mmin, mmax),
            linewidth = 1)
        sinkmap = cm.get_cmap('nickmapSink')
        sinkcolors = sinkmap((np.log10(y) - mmin) / (mmax - mmin))     
        ax.scatter(x, y,marker='.',s=3,facecolor=sinkcolors, lw=0)     
  
  # plot the total mass in sinks too
    if nplotted > 1:
        sel = (np.array(times) <= t)
        x = times[sel]
        y = sinkmasses[sel]
        ax.plot(x, y, color='1.0', linewidth=1.25, alpha=0.6, zorder=1)
  
    ax.set_yscale('log')
    #tmin = np.min(times)
    tmax = np.max(times)
    #tspan = tmax - tmin
    #ax.set_xlim(tmin - 0.1 * tspan, tmax)
    ax.set_xlim(tlimmin, tmax)
    ax.set_ylim(1.0,mlimmax)
    set_ticks(ax, '0.15')
    ax.xaxis.grid(False,which='minor')
    ax.yaxis.grid(False,which='minor')
    
    ax.set_xlabel('time / Myr', fontproperties = tfm, size = 15, color='0.6')
    ax.set_ylabel('sink mass / '+r'M${_\odot}$', fontproperties = tfm, size = 15, color='0.6')
    ax.set_xlabel('time / Myr', fontproperties = tfm, size = 15, color='0.6')
    ax.set_ylabel('sink mass / '+r'M${_\odot}$', fontproperties = tfm, size = 15, color='0.6')  


    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)
    
    plt.savefig(framesdir+'sinkmasses_'+str(snap).zfill(5)+'.png', dpi=400)  
    plt.close()   
                     
                     
