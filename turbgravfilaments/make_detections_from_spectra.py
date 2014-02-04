from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
import shutil
import glob
from astropy.io import ascii
from os.path import expanduser

"""
usage:
run this in the posvel directory that contains all of the spectra.hdf5 files.

this goes through each spectra file and reduced the spectra down to a list of 
'detected' velocities.

e.g. in posvel_0:
python make_detections_from_spectra.py
"""

res = 1024
dres = 1 / res

f = h5py.File('spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

# make empty array of detections
dets = np.ndarray([res, res, 32])
dets.fill(np.nan)

for i in xrange(231,res):
    print i
    ypt = (i + 0.5) * dres
    # read in this row's hdf5 spectra file
    specfile = 'spectra_'+str(i).zfill(4)+'.hdf5'
    f = h5py.File(specfile)
    specs = f['spectra']
    
    # check to make sure we are where we think we are
    if ypt != specs.attrs['z']:
        print 'PROBLEM ',i, ypt, specs.attrs['z']
    
    for s in xrange(res):
        spec = np.array(specs[s])
        # remove anything that was recorded by density = 1.e-99 material
        spec[spec < 1.e-90] = 0
        # normalize the spectrum by its rms
        if sum(spec) > 0:
            spec /= np.sqrt(np.mean(spec**2))
            # set non-'detected' vels equal to zero
            spec[spec < 2] = 0
        
        vcount = 0
        # march through and merge together touching detections
        runningn = 0
        runningv = 0
        if spec[0] >= 1:
            runningn = spec[0]
            runningv = vels[0] * spec[0]
        for v in xrange(1, len(spec)):
            if (spec[v] > 0) & (spec[v-1] == 0):
                # we are starting a new detection
                runningn = spec[v]
                runningv = vels[v] * spec[v]
            if (spec[v] > 0) & (spec[v-1] > 0):
                # we are continuing a detection
                runningn += spec[v]
                runningv += vels[v] * spec[v]
            if (spec[v] == 0) & (spec[v-1] > 0):
                # we are terminating a detection, log it
                dets[i,s,vcount] = runningv / runningn
                vcount += 1
            if vcount > 10:
                print s
                print dets[i,s]
    f.close()
    
    
# NOTE this is [z, y, [vels]]
f = h5py.File('detections.hdf5', 'w')
dset = f.create_dataset('veldetections', data = dets)
f.close()
    