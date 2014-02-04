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
from astropy.io import ascii



# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# the limits on the surface density color map is set up for the compact clouds. here we
# see if we're looking at diffuse clouds, and if so we adjust for that.
snapstr = str(int(sys.argv[1])).zfill(5)
infoname = 'reduced_'+snapstr+'/info_'+snapstr+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)
if boxlen > 5:
    sdoff = np.log10(4)
    vdoff = np.log10(8)
else:
    sdoff = 0.0
    vdoff = 0.0


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    # first do projections
    imshowmap = 'nickmapSD'
    #imshowmap = 'bone_r'
    cdmin = -4.11 - sdoff
    cdmax = 0.069 - sdoff
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.csv'
    
    for i in xrange(3):
        file = fileprefix+'surface_density_'+str(i)+'.hdf5'
        if os.path.exists(file):
            print snap,file
            f = h5py.File(file, 'r')
            sd = f['surface_density']
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
        
            framesdir = outdir+'surfacedensity'+str(i)+'/'
            if not os.path.exists(framesdir):
                os.makedirs(framesdir)
        
            framename = framesdir+'sd'+str(i)+'_frame_'+str(snap).zfill(5)+'.png'
            plt.savefig(framename, dpi = 200)
            f.close() 
            plt.close() 
            del(f)
            del(sd)
            gc.collect()
        
    # now do slices
    imshowmap = 'bone_r'
    imshowmap = 'nickmapVD'
    cdmin = -24.5 - vdoff
    cdmax = -19. - vdoff
    
    for i in xrange(3):
        for j in xrange(5):
            file = fileprefix+'density_slice_'+str(i)+'_'+str(j)+'.hdf5'
            if os.path.exists(file):
                print snap,file
                f = h5py.File(file, 'r')
                sd = f['volume_density']
                print np.min(sd),np.max(sd)
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
        
                framesdir = outdir+'slice'+str(i)+'_'+str(j)+'/'
                if not os.path.exists(framesdir):
                    os.makedirs(framesdir)
        
                framename = framesdir+'sl'+str(i)+'_'+str(j)+'_frame_'+str(snap).zfill(5)+'.png'
                plt.savefig(framename, dpi = 200)
                f.close() 
                plt.close() 
                del(f)
                del(sd)
                gc.collect()    
 

                     
                     
                     
                     
