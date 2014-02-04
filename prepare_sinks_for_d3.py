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
from astropy.table import Table


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')

nsinkmax = 1000
sinktimes = np.zeros([nsinkmax, int(sys.argv[2]) - int(sys.argv[1])])
sinkmasses = np.zeros([nsinkmax, int(sys.argv[2]) - int(sys.argv[1])])

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.csv'
    
    (time, unit_t) = get_time(infoname)
    # convert to Myr
    tmyr = time * unit_t / (31557600 * 1.e6)
    print 'time = ',tmyr
    # see if we have any sink particles to plot
    try:
        with open(sinkname): 
            sinks = ascii.read(sinkname, names=sinkcolumnnames, converters=sinkconverters)
            sinkid = sinks['ID']
            sinkmass = sinks['mass']
            sinktimes[:, snap - int(sys.argv[1])] = tmyr
            for i in xrange(len(sinkid)):
                sinkmasses[sinkid[i], snap - int(sys.argv[1])] = sinkmass[i]
    except IOError:
        pass     
        
sinkdir = './sinkmasstime/'        
if not os.path.exists(sinkdir):
    os.makedirs(sinkdir)        

f = open(sinkdir+'sinkmasses.json', 'w')
f.write('[\n')
# write the first sink
if sinkmasses[1].sum() > 0:
    f.write('{\n')
    f.write('"name": "'+str(1)+'",\n')
    f.write('"data":[ ')
    for ii in xrange(0, len(sinkmasses[1])): # trim all leading zeros but one
	print sinkmasses[1,ii]
        if sinkmasses[1, ii] > 0: 
            break
    print "for sink 1 start at ",ii-1," that's time ",sinktimes[1,ii-1]
    f.write('{"time":'+str(sinktimes[1,ii-1])+', "mass":'+str(sinkmasses[1,ii-1])+'}')
    for j in xrange(ii,len(sinkmasses[1])):
        if sinkmasses[1, j] == 0: # trim trailing zeros in case of merger
            break
        f.write(',{"time":'+str(sinktimes[1,j])+', "mass":'+str(sinkmasses[1,j])+'}')
    f.write(' ]\n')
    f.write('}')
for i in xrange(2,nsinkmax):
    if sinkmasses[i].sum() > 0:
        f.write(', {\n')
        f.write('"name": "'+str(i)+'",\n')
        f.write('"data":[ ')
        for ii in xrange(0, len(sinkmasses[i])): # trim all leading zeros but one
            if sinkmasses[i, ii] > 0: 
                break
        print "for sink ",i," start at ",ii-1," that's time ",sinktimes[i,ii-1]
        f.write('{"time":'+str(sinktimes[i,ii-1])+', "mass":'+str(sinkmasses[i,ii-1])+'}')
        for j in xrange(ii,len(sinkmasses[i])):
            if sinkmasses[i, j] == 0: # trim trailing zeros in case of merger
                break
            f.write(',{"time":'+str(sinktimes[i,j])+', "mass":'+str(sinkmasses[i,j])+'}')
        f.write(' ]\n')
        f.write('}')
f.write(']\n')
f.close()

for i in xrange(nsinkmax):
    if sinkmasses[i].sum() > 0:
        sinkfile = sinkdir + 'sink'+str(i)+'.csv'
        print sinkfile
        data = Table([sinktimes[i], sinkmasses[i]], names=['t', 'm'])
        ascii.write(data, sinkfile, delimiter=',')
        
