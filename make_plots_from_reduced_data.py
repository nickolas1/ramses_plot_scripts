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
    sdoff = np.log10(4)
    vdoff = np.log10(8)
else:
    sdoff = 0.0
    vdoff = 0.0


def make_my_cmaps():
    # surface density color map
    x = [-4.11 - sdoff, -3.11 - sdoff, -.931 - sdoff, .069 - sdoff]
    # x[3] and x[0] are cdmax and cdmin below.
    beginx = (x[1] - x[0]) / (x[3] - x[0])
    begingray = 0.9
    transitionx = (x[2] - x[0]) / (x[3] - x[0])
    transitiongray = 0.35
    finishr = 37 / 256
    finishg = 49 / 256
    finishb = 111 / 256
    cdict = {'red': ((0.0, 1.0, 1.0),
                     (beginx, begingray, begingray),
                     (transitionx, transitiongray, transitiongray),
                     (1.0, finishr, finishr)),
           'green': ((0.0, 1.0, 1.0),
                     (beginx, begingray, begingray),
                     (transitionx, transitiongray, transitiongray),
                     (1.0, finishg, finishg)),
            'blue': ((0.0, 1.0, 1.0),
                     (beginx, begingray, begingray),
                     (transitionx, transitiongray, transitiongray),
                     (1.0, finishb, finishb))}      
    cmap1 = col.LinearSegmentedColormap('my_colormapSD', cdict, N=256, gamma=1.0)
    cm.register_cmap(name='nickmapSD', cmap=cmap1)
    
    
make_my_cmaps()


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
This first part makes plots of the surface density PDF
"""

# override the defaults for this movie plot
mpl.rc('grid', color='0.15')
mpl.rc('grid', linewidth='0.8')
mpl.rc('axes', facecolor='0.0')
mpl.rc('xtick', color='0.6')
mpl.rc('ytick', color='0.6')
mpl.rc('figure', facecolor='0.0')
mpl.rc('savefig', facecolor='0.0')

mu = 2.33 # mean molecular weight
mH = 1.6733e-24

snapstart = int(sys.argv[1])
snapend = int(sys.argv[2])
snapiter = int(sys.argv[3])

for snap in range(snapstart, snapend, snapiter):
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    framesdir = outdir+'surfacedensity0pdfs/'
    if not os.path.exists(framesdir):
            os.makedirs(framesdir)
            
    fig = plt.figure(figsize = (5, 3.5))
    ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
    
    ax.set_xlim(0.99e19,1.01e24)
    ax.set_ylim(0.99e-5,1.01e-1)
    ax.set_yscale('log',nonposy='clip')
    ax.set_xscale('log')
    ax.set_axisbelow(False)
    
    cdensmax = np.log10(10**2.0 / (mu * mH))  # convert these from g to n
    cdensmin = np.log10(10**-6.0 / (mu * mH))
    bins = 128
    binvals = np.arange(cdensmin, 1.000001*cdensmax, (cdensmax - cdensmin) / (bins))
    binmids = 0.5 * (np.roll(binvals, -1) + binvals)
    binmids = binmids[:len(binmids) - 1]
    
    files = [
        'surface_density_0.hdf5',
        'surface_density_1.hdf5',
        'surface_density_2.hdf5']
    files = ['surface_density_0.hdf5']
    colors = [c1,c2,c3]
    
    for i in xrange(len(files)):
        f = h5py.File('reduced_'+str(snap).zfill(5)+'/'+files[i], 'r')
        sd = f['surface_density']
        totalhist = np.zeros(bins)
        print snap
        for j in xrange(sd.shape[0]):
            coldensvals = sd[j]
            coldensvals -= np.log10(mu * mH)
            #if j == 600:
            #    print coldensvals
            #    print cdensmin,cdensmax
            hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
            totalhist += hist  
        f.close() 
        ax.plot(10**binmids, totalhist/np.sum(totalhist), color = '1.0', linewidth = 1.5, alpha=0.8)
    
    plotlim = plt.xlim() + plt.ylim()
    ax.imshow([[0,0],[1,1]], cmap='bone_r', interpolation='bicubic', extent=plotlim,alpha=0.2,zorder=-1)
    ax.fill_between(10**binmids, totalhist/np.sum(totalhist), 10, facecolor='0.0')

    set_ticks(ax, '0.15')
    ax.xaxis.grid(False,which='minor')
    ax.yaxis.grid(False,which='minor')

    ax.set_xlabel(r'surface density / $\mathdefault{cm^{-2}}$', fontproperties = tfm, size = 15, color='0.6')
    #ax.set_ylabel('d', fontproperties = tfm, size = 15)
    ax.set_ylabel(r'volume weighted PDF', fontproperties = tfm, size = 15, color='0.6')

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)
        
    [line.set_zorder(30) for line in ax.lines]
    [line.set_zorder(30) for line in ax.lines]
    
    plt.savefig(framesdir+'SurfaceDensityPDF_'+str(snap).zfill(5)+'.png', dpi=400)   
    plt.close() 


  
                     
"""
This second part makes plots of the sink particle masses with time
"""

# override the defaults for this movie plot
mpl.rc('grid', color='0.15')
mpl.rc('grid', linewidth='0.8')
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


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    
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
    
# plot individual sink masses
mmin = np.log10(1)
mmax = np.log10(100)

for snap in range(snapstart, snapend, snapiter):
    fig = plt.figure(figsize = (5, 3.5))
    ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
    t = alltimes[snap-1]
    print 'plotting sink masses for time ',t
    for i in xrange(int(max(nmt[0,:]))):
        thisone = nmt[:,nmt[0,:] == i]
        if len(thisone[0]) == 0: continue
        x = thisone[2,:]
        if np.min(x) > t: continue
        y = thisone[1,:]
   # if len(x) < 50 and len(x) > 1:
   #     xnew = np.linspace(np.min(x), np.max(x), 50)
   #     ynew = np.interp(xnew, x, y)       
   #     x = xnew
   #     y = ynew
        colorline(x, y, 
            np.log10(y), 
            color='0.8', 
            norm=plt.Normalize(mmin, mmax),
            linewidth = 1)
        sinkmap = cm.get_cmap('nickmapSink')
        sinkcolors = sinkmap((np.log10(y) - mmin) / (mmax - mmin))     
        ax.scatter(x, y,marker='.',s=3,facecolor=sinkcolors, lw=0)     

    ax.set_yscale('log')
    ax.set_xlim(0,2.0)
    ax.set_ylim(0.1,250)
    set_ticks(ax, '0.15')
    ax.xaxis.grid(False,which='minor')
    ax.yaxis.grid(False,which='minor')

    ax.set_xlabel('time / Myr', fontproperties = tfm, size = 15, color='0.6')
    ax.set_ylabel('total sink mass / '+r'M${_\odot}$', fontproperties = tfm, size = 15, color='0.6')

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)

    
    plt.savefig(framesdir+'sinkmasses_'+str(snap).zfill(5)+'.png', dpi=400)  
    plt.close()   
                     
                     
