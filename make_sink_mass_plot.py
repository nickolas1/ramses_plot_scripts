from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys

# import ramses helper functions and get figure directory
homedir = '/home/moon/moeckel/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = '/home/moon/moeckel/Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)

fig = plt.figure(figsize = (5,3.5))

ax = fig.add_axes([0.2, 0.15, 0.75, 0.8])
#ax = fig.add_axes([0., 0., 1., 1.])

times = [0.0]
sinkmasses = [0.0]

indivnames = []
indivmasses = []
indivtimes = []

for snap in range(int(sys.argv[1]),int(sys.argv[2]),1):
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    # see if we have any sink particles to plot
    try:
        with open(sinkname): 
            sinks = get_sinks(sinkname)
            if len(sinks) > 0:
                sinkmass = sinks[:,1].sum() # total mass of sinks in Msun
                (time, unit_t) = get_time(infoname)
		timeMyr = time * unit_t / 31557600.0 / 1.e6
		print snap, timeMyr, sinkmass
                times.append(timeMyr)
                sinkmasses.append(sinkmass)
                timerow = np.ones((sinks.shape[0],1)) * timeMyr
                # nmt = name, mass, time
                nmt = np.hstack((sinks[:,[0,1]],timerow))
                indivnames.append(nmt[:,0])
                indivmasses.append(nmt[:,1])
                indivtimes.append(nmt[:,2])
    except IOError:
        pass      

# flatten lists of individual sink properties, stick in an array
indivnames = [j for i in indivnames for j in i]
indivmasses = [j for i in indivmasses for j in i] 
indivtimes = [j for i in indivtimes for j in i]     
nmt = np.array([indivnames, indivmasses, indivtimes])

# plot total sink mass
ax.plot(times, sinkmasses, color = '0.2', linewidth = 1.5)

# plot individual sink masses
for i in xrange(int(max(nmt[0,:]))):
    thisone = nmt[:,nmt[0,:] == i]
    ax.plot(thisone[2,:], thisone[1,:], color = '0.2', linewidth = 1.25, alpha = 0.5)
    
ax.set_yscale('log')
ax.set_xlim(-.05,7.05)
ax.set_ylim(1,5000)
set_ticks(ax, '0.6')

ax.set_xlabel('time / Myr', fontproperties = tfm, size = 15)
ax.set_ylabel('total sink mass / '+r'M${_\odot}$', fontproperties = tfm, size = 15)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)

    
plt.savefig(outdir+'totalsinkmass.pdf')

    
