from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Circle 
import gc
import sys
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = '/home/moon/moeckel/'
homedir = expanduser('~')+'/'
print homedir
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

# choose dark or light background
background = 'light'

if background == 'dark':
    imshowmap = 'gray'
    textcolor = c3l
    textalpha = 0.5
if background == 'light':
    imshowmap = 'gray_r'
    textcolor = textlightbg
    textalpha = 0.7


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    framename = outdir+'framesD/frame_'+str(snap).zfill(4)+'.png'
    
    projaxis = 0

    (boxlen, unit_l) = get_boxsize(infoname)
    # set column density limits so that the images appear the same.
    # low dens cloud has boxsize 10, high dens cloud has boxsize 5.
    # density of high dens is factor of 8 higher, pathlength factor of 2 smaller.
    # offsets in log of column density limits are thus log10(4)
    if boxlen > 7:
        cdmin = -3.4#-5.1
        cdmax = -0.8#-2.0
    else:
        cdmin = -2.797 #-4.097
        cdmax = -.193

    pf = load(infoname,fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    # choose what region to project
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 1.0 # as a fraction of the whole domain
    wd = 1.0
    res = (720,1280)
    res = (1080,1920)
    
    (lmin, lmax) = get_level_min_max(infoname)
    
    # get a projection of density along x axis
    proj = pf.h.proj('Density', projaxis)
    
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')
    frb = proj.to_frb(width, res, center = cntr, height = height)
    
    fig = plt.figure(figsize = (res[1]/200, res[0]/200), dpi=200)

    #ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])
    ax = fig.add_axes([0., 0., 1., 1.])
    
    # get log of data
    frb_logrho = np.log10(frb['Density'])
    print np.min(frb_logrho), np.max(frb_logrho) 
    dd = pf.h.all_data()
    print dd.quantities['TotalQuantity']('CellMassMsun')

    ax.imshow(frb_logrho,
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
                mincell = 1.0/2**lmax
                sinkrad = 1.5 * mincell
                
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
                sinkrad *= res[1] / wd
                print sinkpos           
                ax.autoscale(False)
                #for s in xrange(len(sinks)):
                #    ax.add_artist(Circle((sinkpos[s,0],sinkpos[s,1]),sinkrad,fc=csink))
	        ax.scatter(sinkpos[:,0],sinkpos[:,1],marker='.',s=4,facecolor=csink,edgecolor=csink)          
    except IOError:
        pass      
        
    # add a time counter
    (time, unit_t) = get_time(infoname)
    timeMyr = time * unit_t / 31557600.0 / 1.e6
    ax.text(.2, .9, r'%.1f' %timeMyr, transform=ax.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=textcolor, alpha=textalpha, snap = False, size=15)
    ax.text(.21, .9, r'Myr', transform=ax.transAxes,
         va = 'baseline',ha='left',fontproperties = lfm, color=textcolor, alpha=textalpha, snap = False, size=13)
    # add a scale bar
    # bar coordinates in units from 0 to 1
    nscales = 2.5  # length of scale bar in units of 0.1 of domain
    barleft = 0.5
    barright = barleft + nscales * 0.1
    barheight = 0.03
    
    #ax2 = ax.twinx() 
    ax2 = fig.add_axes([0,0,.999,.999])
    ax2.plot([barleft, barright], [barheight, barheight], color=textcolor, alpha=textalpha, lw=1.25)
    ax2.patch.set_facecolor('None')
 #   ax2.patch.set_alpha(.1)
    ax2.set_frame_on(False)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    scalePc = nscales * unit_l / 3.086e18 * (boxlen / 10.0) # one unit of length is 10 pc.
    
    if scalePc % 1 > 0.1:
        ax2.text(barright - 0.01 - 0.005, barheight + 0.02, r'%.1f' %scalePc, transform=ax2.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=textcolor, alpha=textalpha, snap = False, size=15) 
    else:
        ax2.text(barright - 0.01 - 0.005, barheight + 0.02, r'%.0f' %scalePc, transform=ax2.transAxes,
         va = 'baseline',ha='right',fontproperties = lfm, color=textcolor, alpha=textalpha, snap = False, size=15) 
    ax2.text(barright - 0.01 + 0.005, barheight + 0.02, 'pc', transform=ax2.transAxes,
         va = 'baseline',ha='left',fontproperties = lfm, color=textcolor, alpha=textalpha, snap = False, size=13)
    
      
    
        
    plt.savefig(framename, dpi = 200)
    
    del(frb)
    del(frb_logrho)
    del(pf)
    gc.collect()
    
