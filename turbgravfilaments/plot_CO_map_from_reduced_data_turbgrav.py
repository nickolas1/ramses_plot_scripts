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
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from matplotlib import rcParams


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12) 

rcParams['xtick.direction'] = 'out'

outdir = get_output_path(homedir)
outdir = './'
# the limits on the surface density color map is set up for the compact clouds. here we
# see if we're looking at diffuse clouds, and if so we adjust for that.
snapstr = str(int(sys.argv[1])).zfill(5)
infoname = 'reduced_'+snapstr+'/info_'+snapstr+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    # first do projections
    imshowmap = 'nickmapVD2'
    #imshowmap = 'bone_r'
    
    """
    this data is generated by making a surface density map taking into account only
    gas that is between 10^3 and 10^4.5 cm^-3. to convert to a crude approximation 
    of a C18O map, use
    David S. Meier and Jean L. Turner ApJ 551:687 2001 equation 2
    
    N(H2)C18O = 2.42e14 cm^-2 [H2]/[C18O] * exp(5.27/Tex)/(exp(5.27/Tex)-1) IC18O K km/s
    [H2]/[C18O] = 2.94e6
    
    so first convert g cm^-2 to cm^-2 using mu = 2.33, then convert to ICO using Tex=10 K
    """
    cdmin = 10**-3.3 
    cdmax = 10**-1.5
    cdmin = 0
    cdmax = 4
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.csv'
    
    for i in xrange(3):
        file = fileprefix+'surface_density_CO'+str(i)+'.hdf5'
        if os.path.exists(file):
            print snap,file
            f = h5py.File(file, 'r')
            sd = f['surface_density_CO']
            sd = 10**np.array(sd)  # convert to linear units
            print np.mean(sd),np.max(sd)
            sd /= (2.33 * const.m_p.cgs.value) # convert to number density
            print np.mean(sd),np.max(sd)
            sd /= (2.42e14 * 2.94e6) # non-temperature factors of IC18O conversion
            sd /= (np.exp(5.27/10) / (np.exp(5.27/10) - 1)) # temperature part
            print np.mean(sd),np.max(sd)
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
            ax2.set_xlabel(r'$\mathdefault{I_{C^{18}O}}$ / K km s$\mathdefault{^{-1}}$', fontproperties = tfm, size=15, color='0.15')
            
            set_ticks(ax2, '0.15')
            for label in ax2.get_xticklabels() + ax2.get_yticklabels():
                label.set_fontproperties(tfm)


            """
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
               add a scalebar
               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            """
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
                va = 'bottom', ha = 'left', fontproperties = tfm, color='0.15', snap = False)
            set_ticks(ax3, '0.15')
            for label in ax3.get_xticklabels() + ax3.get_yticklabels():
                label.set_fontproperties(tfm)
                
                
                            
            framesdir = outdir+'surfacedensityCO'+str(i)+'/'
            if not os.path.exists(framesdir):
                os.makedirs(framesdir)
        
            framename = framesdir+'sd'+str(i)+'_frame_'+str(snap).zfill(5)+'.png'
            plt.savefig(framename, dpi = 200)
            f.close() 
            plt.close() 
            del(f)
            del(sd)
            gc.collect()

