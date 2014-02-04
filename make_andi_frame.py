from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
from astropy.io import ascii
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font main
    fname=fontdir+'Gotham-Book.ttf', size=13)    
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=11)  


for snap in range(int(sys.argv[1]),int(sys.argv[2]),1):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    framename = outdir+'framesD/frame_'+str(snap).zfill(4)+'.png'
    
    projaxis = 1

    (boxlen, unit_l) = get_boxsize(infoname)
    # set column density limits so that the images appear the same.
    # low dens cloud has boxsize 10, high dens cloud has boxsize 5.
    # offsets in log of column density limits are thus log10(8)
    if boxlen > 7:
        cdmin = -5.1
        cdmax = -2.0
    else:
        cdmin = -4.19
        cdmax = -1.09

    pf = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    # get a projection of density along y axis
    proj = pf.h.proj('Density', 1)

    wd = 0.5 # this is messed up- figure it out. yt might not get size right.
    # res should be base resolution times 2**levels of refinement * wd
    res = (512,512)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    frb = proj.to_frb(width, res, center = cntr, height = height)
    ascii.write(np.log10(frb['Density']),'snap72_surfacedensity_xz.ascii')
    
    
    # get a projection of density along x axis
    proj = pf.h.proj('Density', 0)
    frb = proj.to_frb(width, res, center = cntr, height = height)
    ascii.write(np.log10(frb['Density']),'snap72_surfacedensity_yz.ascii')
    
    
    # get a projection of density along z axis
    proj = pf.h.proj('Density', 2)
    frb = proj.to_frb(width, res, center = cntr, height = height)
    ascii.write(np.log10(frb['Density']),'snap72_surfacedensity_xy.ascii')
    
    
    sys.exit()
    
    fig = plt.figure(figsize = (res[1]/200, res[0]/200), dpi=200)

    #ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])
    ax = fig.add_axes([0., 0., 1., 1.])
    
    # get log of data
    frb_logrho = np.log10(frb['Density'])
    
    ax.imshow(frb_logrho,
              interpolation = 'nearest',
              origin = 'lower',
              vmin = cdmin,
              vmax = cdmax,
              cmap = 'bone')
    
    # see if we have any sink particles to plot
    try:
        with open(sinkname): 
            sinks = get_sinks(sinkname)
            if len(sinks) > 0:
                sinkpos = sinks[:,2:5]
                #sinkpos = sinks[:,1:4]
                sinkpos[:] /= boxlen # shrink to 0-1 in all dimensions
                # get projected positions
                keep = np.array([1,1,1])
                keep[projaxis] = 0
                keep = np.array(keep, dtype=bool)
                sinkpos = sinkpos[:,keep]
            
                # restrict to same region as density plot
                ledge = cntr[0] - wd/2
                bedge = cntr[2] - ht/2
                sinkpos[:] -= np.array([ledge, bedge])
                # convert to imshow scale
                sinkpos *= res[1] / wd
                print sinkpos           
                ax.autoscale(False)
                ax.scatter(sinkpos[:,0],sinkpos[:,1],marker='.',s=6,facecolor='r',edgecolor='r')          
    except IOError:
        pass      
        
    # add a time counter
    (time, unit_t) = get_time(infoname)
    timeMyr = time * unit_t / 31557600.0 / 1.e6
    ax.text(.2, .9, r'%.1f' %timeMyr, transform=ax.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=c3l, alpha=0.5, snap = False, size=15)
    ax.text(.21, .9, r'Myr', transform=ax.transAxes,
         va = 'baseline',ha='left',fontproperties = lfm, color=c3l, alpha=0.5, snap = False, size=13)
    # add a scale bar
    # bar coordinates in units from 0 to 1
    nscales = 2.5  # length of scale bar in units of 0.1 of domain
    barleft = 0.5
    barright = barleft + nscales * 0.1
    barheight = 0.07
    
    #ax2 = ax.twinx() 
    ax2 = fig.add_axes([0,0,.999,.999])
    ax2.plot([barleft, barright], [barheight, barheight], color=c3l, alpha=0.5, lw=1.25)
    ax2.patch.set_facecolor('None')
 #   ax2.patch.set_alpha(.1)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    scalePc = nscales * unit_l / 3.086e18 * (boxlen / 10.0) # one unit of length is 10 pc.
    
    if scalePc % 1 > 0.1:
        ax2.text(barright - 0.01 - 0.005, barheight + 0.02, r'%.1f' %scalePc, transform=ax2.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=c3l, alpha=0.5, snap = False, size=15) 
    else:
        ax2.text(barright - 0.01 - 0.005, barheight + 0.02, r'%.0f' %scalePc, transform=ax2.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=c3l, alpha=0.5, snap = False, size=15) 
    ax2.text(barright - 0.01 + 0.005, barheight + 0.02, 'pc', transform=ax2.transAxes,
         va = 'baseline',ha='left',fontproperties = lfm, color=c3l, alpha=0.5, snap = False, size=13)
    
      
    
        
    plt.savefig(framename, dpi = 200)
    
    del(frb)
    del(frb_logrho)
    del(pf)
    gc.collect()
    
