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


def make_my_cmap():
    x = [-4.11, -3.11, -.931, .069]
    # x[3] and x[0] are cdmax and cdmin below.
    beginx = (x[1] - x[0]) / (x[3] - x[0])
    begingray = 0.9
    transitionx = (x[2] - x[0]) / (x[3] - x[0])
    transitiongray = 0.35
    finishr = 37.0 / 256
    finishg = 49.0 / 256
    finishb = 111. / 256
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
    cmap1 = col.LinearSegmentedColormap('my_colormap', cdict, N=256, gamma=1.0)
    cm.register_cmap(name='nickmap', cmap=cmap1)

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
    
make_my_cmap()

# convert these to 
mu = 2.33 # mean molecular weight
mH = 1.6733e-24

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    fig = plt.figure(figsize = (5, 3.5))
    ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
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
    colors = [c1,c2,c3]
    
    imshowmap = 'nickmap'
    cdmin = -4.11 
    cdmax = 0.069
    
    print snap
    
    for i in xrange(len(files)):
        f = h5py.File('output_'+str(snap).zfill(5)+'/'+files[i], 'r')
        sd = f['surface_density']
        print np.max(sd)
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
        
        framename = outdir+'frames'+str(i).zfill(1)+'/frame_'+str(snap).zfill(4)+'.png'
        
        plt.savefig(framename, dpi = 200)
        del(fig)
        gc.collect()
 

                     
                     
                     
                     