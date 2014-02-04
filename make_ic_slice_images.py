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


for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    fileprefix = 'ic_slices_reduced_'+str(snap).zfill(5)+'/'
    
    framesdir = outdir+'ic_slices/'
    if not os.path.exists(framesdir):
        os.makedirs(framesdir)
    
    i = 0
    
    # now do slices
    imshowmap = 'bone_r'
    imshowmap = 'nickmapVD'
    vdoff = np.log10(8)
    cdmin = -24.5 - vdoff 
    cdmax = -19. - vdoff 
   # imshowmap = 'gray_r'
    cdmin = - 24
    cdmax = -20.25
    
    for j in xrange(256,512):
        file = fileprefix+'density_slice_'+str(i)+'_'+str(j).zfill(5)+'.hdf5'
        print snap,file
        f = h5py.File(file, 'r')
        sd = f['volume_density']
        print np.min(sd),np.max(sd), cdmin, cdmax
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

        framename = framesdir+'sl'+str(i)+'_'+str(j)+'_frame_'+str(snap).zfill(5)+'.png'
        plt.savefig(framename, dpi = 200)
        f.close() 
        plt.close() 
        del(f)
        del(sd)
        gc.collect()    
 

                     
                     
                     
                     
