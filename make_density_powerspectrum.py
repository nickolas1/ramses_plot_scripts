from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys

# import ramses helper functions and get figure directory
homedir = '/home/moon/moeckel/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font main
    fname=fontdir+'Gotham-Book.ttf', size=13)    
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=11)  


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    framename = outdir+'densityPDF'+str(snap).zfill(4)+'.png'
    
    (boxlen, unit_l) = get_boxsize(infoname)

    ds = load(infoname,fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]
    
    sphererad = 27.0
    spherevol = 4.0 * np.pi / 3.0 * (sphererad * 3.086e18)**3
    sp = ds.h.sphere(cntr, (sphererad, 'pc'))
    
    gasmass = sp.quantities['TotalQuantity']('CellMassMsun')
    print gasmass
    gasmass = sp.quantities['TotalQuantity']('CellMassCode')
    print gasmass
    
    nbins = 128
    dmin = 1.e-25
    dmax = 1.e-18

    profile = BinnedProfile1D(sp,nbins,'Density',dmin,dmax)
    profile.add_fields("CellVolume", weight=None)
    
    fig = plt.figure(figsize = (5,3.5))

    ax = fig.add_axes([0.2, 0.15, 0.75, 0.8])
   
    # plot density pdf
    ax.plot(profile['Density'], profile['CellVolume']/spherevol, color = '0.2', linewidth = 1.5)
    ax.yaxis.grid(False,which='minor')   
    ax.xaxis.grid(False,which='minor')
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_xlim(-.05,7.05)
    ax.set_ylim(1.e-8,0.1)
    set_ticks(ax, '0.6')
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)
    
    plt.savefig(framename, dpi = 200)
    
    
