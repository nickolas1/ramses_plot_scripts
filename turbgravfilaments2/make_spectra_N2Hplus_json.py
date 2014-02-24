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
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from scipy import special
from mpl_toolkits.mplot3d import Axes3D

"""
usage:
python fooppv.py N A F
N: number of output to use. reduced_N needs to be here.
A: axis of the projection (0, 1, 2)
"""

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *


# downsample the spectra to a 256 squared grid
inres = 1024#1024
outres = 256#256
stride = int(inres / outres)
print stride

snap = int(sys.argv[1])
axis = int(sys.argv[2])
fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

f = h5py.File(fileprefix+'posvel_'+str(axis)+'/spectrumvels.hdf5')
vels = np.array(f['binmidskms'])
f.close()

f = open('velocities.json', 'w')
f.write('[\n{"velocity":[')
for v in xrange(len(vels)):
    f.write(str(round(vels[v],3)))
    if v < len(vels) - 1: 
        f.write(',')
f.write(']}\n]')
f.close()

outspecs = np.array([[np.zeros(len(vels)) for x in xrange(outres)] for y in xrange(outres)])

#
# j is image up
# |
# | 
# |_______i is image right
#
# spectra files march along j
#
totnonzero = 0
for inj in xrange(inres):
    specfile = fileprefix+'posvel_'+str(axis)+'/spectra_N2Hplus_'+str(inj).zfill(4)+'.hdf5'
    f = h5py.File(specfile)
    specs = f['spectraN2Hplus']
    print specfile
    outj = inj//stride
    for ini in xrange(inres):
        outi = ini//stride
        outspecs[outi,outj] += specs[ini]
        totnonzero += len(specs[ini][specs[ini]>0])
    f.close()
meannonzero = outspecs.sum() / totnonzero
#outspecs /= meannonzero
nonzero = outspecs[outspecs > 0]

print np.mean(outspecs),' ',np.median(outspecs),' ',np.max(outspecs)
print np.mean(nonzero),' ',np.median(nonzero),' ',np.max(nonzero)

"""
outspecs /= np.median(nonzero)
print np.mean(outspecs),' ',np.median(outspecs),' ',np.max(outspecs)

# if doing this, round to 1 decimal in the write section below

# median is 3.15197674139e-24 for 256
# median is 1.0888351136e-23 for 128. use the same conversion for n2h+?
"""

# new strategy: make the highest value = 512, divide by that and round.
# highest density for snapshot_00018 in C18O is 6.43825472928e-21
outspecs /= 6.43825472928e-21
outspecs *= 512
# if doing this, int the rounded values in the write section below
outspecs = np.round(outspecs,0)

f = open('spectra_N2Hplus_'+str(snap).zfill(5)+'.json', 'w')
f.write('[\n')
for j in xrange(outres):
    for i in xrange(outres): 
        if outspecs[i,j].sum() > 0.2:
            f.write('{"c":'+str(i + j*outres)+',"s":[')
            for v in xrange(len(vels)):
                if outspecs[i,j,v] > 0:
                    sstr = str(int(outspecs[i,j,v]))
                else:
                    sstr = '0'
                f.write(sstr)
                if v < len(vels) - 1: 
                    f.write(',')      
            f.write(']}')
            if i + j*outres < outres**2 - 1:
                f.write(',\n')
f.write(']')
f.close()

f = open('spectra_N2Hplus_'+str(snap).zfill(5)+'_sparse.json', 'w')
f.write('[\n')
for j in xrange(outres):
    for i in xrange(outres): 
        if outspecs[i,j].sum() > 0.2:
            f.write('{"c":'+str(i + j*outres)+',"s":[')
            zerocount = 0
            for v in xrange(len(vels)):
                val = int(outspecs[i,j,v])
                if val > 0:
                    if zerocount > 0:
                        f.write(str(-zerocount))
                        f.write(',')
                        zerocount = 0
                    f.write(str(val))
                    if v < len(vels) - 1: 
                        f.write(',') 
                else:
                    zerocount += 1    
            if zerocount > 0:
                f.write(str(-zerocount))
            f.write(']}')
            if i + j*outres < outres**2 - 1:
                f.write(',\n')
f.write(']')
f.close()

    

   