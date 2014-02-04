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
if boxlen > 7:
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
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.out'
    
    for i in xrange(3):
        file = fileprefix+'surface_density_CO_'+str(i)+'.hdf5'
        if os.path.exists(file):
            print snap,file
            f = h5py.File(file, 'r')
            sd = f['surface_density_CO']
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
                    sinks = get_sinks(sinkname)
                    if len(sinks) > 0:
                        # figure out the size of the sinks in units of 0-1
                        #mincell = 1.0/2**lmax
                        #sinkrad = 1.5 * mincell
                
                        sinkpos = sinks[:,2:5]
                        sinkpos[:] /= boxlen # shrink to 0-1 in all dimensions
                        # get projected positions
                        keep = np.array([1,1,1])
                        keep[i] = 0
                        keep = np.array(keep, dtype=bool)
                        sinkpos = sinkpos[:,keep]
            
                        # restrict to same region as density plot
                        #ledge = cntr[0] - wd/2
                        #bedge = cntr[2] - ht/2
                        #sinkpos[:] -= np.array([ledge, bedge])
                        # convert to imshow scale
                        #sinkpos *= res[1] / wd
                        #sinkrad *= res[1] / wd
                        sinkpos *= sd.shape[0]
                        print sinkpos      
                        sinkmass = sinks[:,1]
                        # color by the log of mass. the minimum that we plot is 0.1 Msun,
                        # max is a few hundred.
                        mmin = np.log10(1)
                        mmax = np.log10(100)
                        sinkmap = cm.get_cmap('nickmapSink')
                        sinkcolors = sinkmap((np.log10(sinkmass) - mmin) / (mmax - mmin))     
                        ax.autoscale(False)
                        #for s in xrange(len(sinks)):
                        #    ax.add_artist(Circle((sinkpos[s,0],sinkpos[s,1]),sinkrad,fc=csink))
                        ax.scatter(sinkpos[:,0],sinkpos[:,1],marker='.',s=9,facecolor=sinkcolors,lw=0)              
            except IOError:
                pass     
        
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

