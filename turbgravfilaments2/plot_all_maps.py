from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.cm as cm
import gc
import sys
import h5py
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from matplotlib import rcParams
import scipy.ndimage as ndimage


"""
plot column density from different density ranges. 

the hdf5 files contain log surface density in g cm^-2
files with 'C18O': surface density of gas with number density 10^3 < n < 10^4.5
files with 'N2Hplus': surface density of gas with number density 10^4.5 < n
files that just say 'surface_density': all the gas

this script plots column density in cm^-2, converting from g cm^-2 using mu=2.33
"""


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *

#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=7)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=6) 

rcParams['xtick.direction'] = 'out'

outdir = get_output_path(homedir)
#outdir = './'

snapstr = str(int(sys.argv[1])).zfill(5)
infoname = 'reduced_'+snapstr+'/info_'+snapstr+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)


"""
   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   toggle these to plot a colorbar and length scale bar
   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
"""
PlotColorBars = False
PlotScaleBar = False
PlotSinks = False

outputres = 768
outputres = 1024

mu = 2.33

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    # first do projections
    imshowmap = 'nickmapVD2'
    #imshowmap = 'bone_r'
    
    cdmin = 10**-3.3 
    cdmax = 10**-1.5
    cdmin = 0
    cdmax = 4
    
    cdmin = 0
    cdmax = 1.*10**22
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.csv'

    
    for i in xrange(3):
        fileCO = fileprefix+'surface_density_C18O'+str(i)+'.hdf5'
        fileN2Hplus = fileprefix+'surface_density_N2Hplus'+str(i)+'.hdf5'
        fileSD = fileprefix+'surface_density'+str(i)+'.hdf5'
        if os.path.exists(fileCO):
            print ''
            print ''
            print snap,fileCO
            fC18O = h5py.File(fileCO, 'r')
            sdC18O = fC18O['surface_density_C18O']
            # convert to linear units, divide by mu * mH
            sdC18O = 10**np.array(sdC18O) / (mu * const.m_p.cgs.value)
            # there are a lot of very small values- cull them before getting some
            # info on the range of interesting values  
            sdnonzero = sdC18O[sdC18O > 10**4]
            print 'mean non-zero C18O column density: ',np.mean(sdnonzero),'cm^-2'
            print 'median non-zero C18O column density: ',np.median(sdnonzero),'cm^-2'
            print 'max C18O column density: ',np.max(sdC18O),'cm^-2'
            
 
            fN2Hplus = h5py.File(fileN2Hplus, 'r')
            sdN2Hplus = fN2Hplus['surface_density_N2Hplus']
            # convert to linear units, divide by mu * mH
            sdN2Hplus = 10**np.array(sdN2Hplus) / (mu * const.m_p.cgs.value)
            # there are a lot of very small values- cull them before getting some
            # info on the range of interesting values  
            sdnonzero = sdN2Hplus[sdN2Hplus > 10**4]
            print ''
            print 'mean non-zero N2Hplus column density: ',np.mean(sdnonzero),'cm^-2'
            print 'median non-zero N2Hplus column density: ',np.median(sdnonzero),'cm^-2'
            print 'max N2Hplus column density: ',np.max(sdN2Hplus),'cm^-2'
            
            
            fSD = h5py.File(fileSD, 'r')
            sdSD = fSD['surface_density']
            # convert to linear units, divide by mu * mH
            sdSD = 10**np.array(sdSD) / (mu * const.m_p.cgs.value)
            # back to log units
            sdSD = np.log10(sdSD)
            print ''
            print 'min log total column density: ',np.min(sdSD),'cm^-2'
            print 'mean log total column density: ',np.mean(sdSD),'cm^-2'
            print 'median log total column density: ',np.median(sdSD),'cm^-2'
            print 'max log total column density: ',np.max(sdSD),'cm^-2'
            
            
            cdmin = 0
            cdmax = 1.0725e22  # set this to max / 2
            imshowmap = 'nickmapVD2'
            
            fig = plt.figure(figsize = (outputres/128, outputres/128), dpi=128)
            ax = fig.add_axes([0., 0., 1., 1.])
            #sdC18O = np.ma.masked_where(sdC18O < 10**20, sdC18O)
            ax.imshow(sdC18O,
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
            if PlotSinks:
                try:
                    with open(sinkname): 
                        sinks = ascii.read(sinkname, names=sinkcolumnnames, converters=sinkconverters, data_start=0)
                        if len(sinks['ID']) > 0: 
                            # figure out the size of the sinks in units of 0-1
                            #mincell = 1.0/2**lmax
                            #sinkrad = 1.5 * mincell
                            sscale = sd.shape[0] / boxlen
                        
                            if i == 0:
                                i0 = 'y'
                                i1 = 'z'
                            if i == 1:
                                i0 = 'x'
                                i1 = 'z'
                            if i == 2:
                                i0 = 'x'
                                i1 = 'y'

                            # convert to imshow scale
                            #sinkpos *= res[1] / wd
                            #sinkrad *= res[1] / wd

                            # color by the log of mass. the minimum that we plot is 0.1 Msun,
                            # max is a few hundred.
                            mmin = np.log10(1)
                            mmax = np.log10(100)
                            sinkmap = cm.get_cmap('nickmapSink')
                            sinkcolors = sinkmap((np.log10(sinks['mass']) - mmin) / (mmax - mmin))     
                            ax.autoscale(False)
                            #for s in xrange(len(sinks)):
                            #    ax.add_artist(Circle((sinkpos[s,0],sinkpos[s,1]),sinkrad,fc=csink))
                            ax.scatter(sinks[i0]*sscale,sinks[i1]*sscale,marker='.',s=9,facecolor=sinkcolors,lw=0)                     
                except IOError:
                    pass     
                
            """
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
               add a colorbar
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            """
            if PlotColorBars:
                ax2 = fig.add_axes([0.1, 0.1, 0.4, 0.015])
                a = np.outer(np.arange(cdmin, cdmax, (cdmax - cdmin)/255), np.ones(10)).T
                ax2.imshow(a, 
                    aspect = 'auto',
                    interpolation = 'nearest',
                    origin = 'lower',
                    vmin = cdmin,
                    vmax = cdmax,
                    cmap = imshowmap,
                    extent = [cdmin, cdmax, 0, 1])
                ax2.set_frame_on(False)
                ax2.axes.get_yaxis().set_visible(False)
                ax2.xaxis.set_ticks(np.arange(cdmin, cdmax+1, 1.0))
                ax2.set_xlabel(r'$\mathdefault{I_{C^{18}O}}$ / K km s$\mathdefault{^{-1}}$', fontproperties = tfm, size=8, color='0.15')
            
                set_ticks(ax2, '0.15')
                for label in ax2.get_xticklabels() + ax2.get_yticklabels():
                    label.set_fontproperties(tfm)


            """
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
               add a scalebar
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            """
            if PlotScaleBar:
                ax3 = fig.add_axes([0.7, 0.1, 0.2, 0.0015])
                a = np.outer(np.ones(100)*.8*cdmax, np.ones(10)).T
                ax3.imshow(a, 
                    aspect = 'auto',
                    interpolation = 'nearest',
                    origin = 'lower',
                    vmin = cdmin,
                    vmax = cdmax,
                    cmap = imshowmap,
                    extent = [cdmin, cdmax, 0, 1])
                ax3.set_frame_on(False)
                ax3.axes.get_yaxis().set_visible(False)
                ax3.axes.get_xaxis().set_visible(False)
                ax3.text(0.1, 0.75, r'2pc', transform = ax3.transAxes,
                    va = 'bottom', ha = 'left', fontproperties = tfm, size=8, color='0.15', snap = False)
                set_ticks(ax3, '0.15')
                for label in ax3.get_xticklabels() + ax3.get_yticklabels():
                    label.set_fontproperties(tfm)
                
                
                            
            framesdir = outdir+'surfacedensityCO'+str(i)+'/'
            framesdir = './'
            if not os.path.exists(framesdir):
                os.makedirs(framesdir)
            framename = framesdir+'C18O_'+str(i)+'_frame_'+str(snap).zfill(5)+'.png'
            plt.savefig(framename, dpi = 128, transparent = True)
            plt.close()
            
            """
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
               contours of NH2Plus
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            """
        
            fig = plt.figure(figsize = (2*outputres/128, 2*outputres/128), dpi=128)
            ax = fig.add_axes([0., 0., 1., 1.])
            contourcolors = [(252/255, 182/255, 252/255),
                        (247/255, 141/255, 255/255),
                        (246/255, 103/255, 255/255),
                        (241/255, 56/255, 255/255),
                        (236/255, 0/255, 255/255)]
            contourvals = [5.e20, 1.e21, 1.5e21, 2.e21, 2.5e21, 2.e21]
            ax.contour(ndimage.gaussian_filter(sdN2Hplus, sigma=0.5, order=0),contourvals,
                origin = 'lower',
                linewidths = 0.75,
                colors = contourcolors)
            
            # turn off axes
            ax.set_frame_on(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_xaxis().set_visible(False) 
            
            framename = framesdir+'N2Hplus_'+str(i)+'_frame_'+str(snap).zfill(5)+'.png'
            plt.savefig(framename, dpi = 128, transparent = True)
            plt.close()
            
            """
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
               projected total surface density
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            """
            
            cdmin = 20.8
            cdmax = 23.4
            imshowmap = 'nickmapSD2'
            
            fig = plt.figure(figsize = (outputres/128, outputres/128), dpi=128)
            ax = fig.add_axes([0., 0., 1., 1.])
            ax.imshow(sdSD,
                  interpolation = 'nearest',
                  origin = 'lower',
                  vmin = cdmin,
                  vmax = cdmax,
                  cmap = imshowmap)
            # turn off axes
            ax.set_frame_on(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_xaxis().set_visible(False)
            
            framename = framesdir+'totalcolumn_'+str(i)+'_frame_'+str(snap).zfill(5)+'.png'
            plt.savefig(framename, dpi = 128, transparent = True)
            
            fC18O.close() 
            fN2Hplus.close()
            fSD.close()
            plt.close() 
            del(fig)
            del(fC18O)
            del(fN2Hplus)
            del(sdC18O)
            del(sdN2Hplus)
            del(sdnonzero)
            gc.collect()

