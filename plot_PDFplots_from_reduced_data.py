from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
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


"""
This first part makes plots of the surface density PDF
"""

# override the defaults for this movie plot
mpl.rc('grid', color='0.15')
mpl.rc('grid', linewidth='1.0')
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
    print plotlim
    ax.imshow([[0,0],[1,1]], cmap='bone_r', interpolation='bicubic', extent=plotlim,alpha=0.4,zorder=-1)
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


  
                     

