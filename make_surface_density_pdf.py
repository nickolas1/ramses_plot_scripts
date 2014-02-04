from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
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

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)
    


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

    for i in xrange(len(files)):
        f = h5py.File('output_'+str(snap).zfill(5)+'/'+files[i], 'r')
        sd = f['surface_density']
        totalhist = np.zeros(bins)
        for j in xrange(sd.shape[0]):
            coldensvals = sd[j]
            coldensvals -= np.log10(mu * mH)
            if j == 600:
                print coldensvals
                print cdensmin,cdensmax
            hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
            totalhist += hist  
        f.close() 
        print totalhist

        ax.plot(10**binmids, totalhist, color = colors[i], linewidth = 1.5, alpha=0.8)
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlim(5.e18,1.e24)
    ax.set_ylim(1,1.e5)

    set_ticks(ax, '0.6')
    ax.xaxis.grid(False,which='minor')
    ax.yaxis.grid(False,which='minor')

    plotlim = mpl.xlim() + mpl.ylim()
    print plotlim
    ax.imshow([0,0],[1,1], cmap=mp.cm.Grays, interpolation='bicubic', extent=plotlim)

    ax.set_xlabel(r'surface density / $\mathdefault{cm^{-2}}$', fontproperties = tfm, size = 15)
    #ax.set_ylabel('d', fontproperties = tfm, size = 15)
    ax.set_ylabel(r'volume weighted PDF', fontproperties = tfm, size = 15)

    (time, unit_t) = get_time(infoname)
    timeMyr = time * unit_t / 31557600.0 / 1.e6
    horiz = 5.e22
    vert = 2.e4
    ax.text(horiz, vert, r'%.1f' %timeMyr, transform = ax.transData, 
        ha = 'right',va = 'baseline', fontproperties = lfm, color = c1, snap = False)
    ax.text(1.1*horiz, vert, r'Myr', transform = ax.transData,
        ha = 'left', va = 'baseline', fontproperties = lfm, color = c1, snap = False)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)

    plt.savefig(outdir+'SurfaceDensities/SurfaceDensityPDF_'+str(snap).zfill(5)+'.png')   
    plt.clf()
    
    res = (1080,1920)
    fig = plt.figure(figsize = (res[1]/200, res[0]/200), dpi=200)

    i = 0
    f = h5py.File('output_'+str(snap).zfill(5)+'/'+files[i], 'r')
    sd = f['surface_density']
    fig = plt.figure(figsize = (sd.shape[0]/200,sd.shape[0]/200), dpi=200)
    ax = fig.add_axes([0.,0.,1.,1.])
    cdmin = -2.797 #-4.097
    cdmax = -.193
    ax.imshow(sd,
          interpolation = 'nearest',
          origin = 'lower',
          vmin = cdmin,
          vmax = cdmax,
          cmap = 'gray_r')
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    framename = outdir+'framesD/framefoo_'+str(snap).zfill(4)+'.png'
    plt.savefig(framename, dpi = 200)

    
