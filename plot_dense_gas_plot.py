from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
import h5py
from astropy.io import ascii
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)
    
fig = plt.figure(figsize = (5,3.5))

ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
#ax = fig.add_axes([0., 0., 1., 1.])

times = []
mass1 = []
mass2 = []
mass3 = []
mass4 = []
sinkmasses = []


# set limits for density that we're interested in
ndense1 = 1.e3
ndense2 = 1.e4
ndense3 = 1.e5
mu = 2.33 # mean molecular weight
mH = 1.6733e-24

mdense1 = ndense1 * mu * mH
mdense2 = ndense2 * mu * mH
mdense3 = ndense3 * mu * mH

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    dataname = 'reduced_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins.dat'
    infoname = 'reduced_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'reduced_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    limitsname = './densegasplotlimits.dat'
    
    f = open(limitsname)
    line = f.readline()
    sl = line.split()
    alltlow = float(sl[0])
    allthi = float(sl[1])
    line = f.readline()
    sl = line.split()
    qtlow = float(sl[0])
    qthi = float(sl[1])
    f.close()
    
    if os.path.isfile(dataname):
        (time, unit_t) = get_time(infoname)
        timeMyr = time * unit_t / 31557600.0 / 1.e6 + 1.e-10
        times.append(timeMyr)
    
        data = ascii.read(dataname)
        dense1sel = (data['Density'] > mdense1)
        dense2sel = (data['Density'] > mdense2)
        dense3sel = (data['Density'] > mdense3)
    
        cumulativemass = np.cumsum(data['CellMassMsun'])
        massabove = cumulativemass.max() - cumulativemass
        logdens = np.log10(data['Density'])
   
        # should interpolate this linearly in log density probably
        md1 = np.sum(data['CellMassMsun'][dense1sel])
        md2 = np.sum(data['CellMassMsun'][dense2sel])
        md3 = np.sum(data['CellMassMsun'][dense3sel])

        md1 = np.interp(np.log10(mdense1), logdens, massabove)
        md2 = np.interp(np.log10(mdense2), logdens, massabove)
        md3 = np.interp(np.log10(mdense3), logdens, massabove)
    
        mass1.append(md1)
        mass2.append(md2)
        mass3.append(md3)
         
        # see if we have any sink particles 
        sinks = get_sinks(sinkname)
        if len(sinks) > 0:
            sinkmass = np.sum(sinks[:,1]) # total mass of sinks in Msun
        else:
            sinkmass = 0.0
        sinkmasses.append(sinkmass)

     
        print snap, timeMyr, md1, md2, md3, sinkmass
        
    # now get surface density 'dense gas' as well. 
    files = [
            'surface_density_0.hdf5',
            'surface_density_1.hdf5',
            'surface_density_2.hdf5']       
    (boxlen, unit_l) = get_boxsize(infoname)
    (levelmin, levelmax) = get_level_min_max(infoname)
    print boxlen, unit_l
    print levelmin, levelmax
    thism4 = []
    for i in xrange(len(files)):
        f = h5py.File('reduced_'+str(snap).zfill(5)+'/'+files[i], 'r')
        sd = f['surface_density']
        md4 = 0.0
        for j in xrange(sd.shape[0]):
            coldensvals = sd[j]
            # these are log10(g / cm^2). convert to msun / pc^2
            coldensvals = 10**coldensvals
            coldensvals *= 4785.63 # that's pc**2 / msun
            sdenscut = 120.0 # surface density cut in msun / pc^2
            # get total mass above the surface density cut
            md4 += np.sum(coldensvals[coldensvals >= sdenscut])
            # scale to the area of each pixel
            # is 10 pc, and we want this in pc.
        pixlenau = boxlen * unit_l / 2**levelmax / 3.086e18
        md4 *= pixlenau**2
        thism4.append(md4)
        f.close() 
    print thism4
    mass4.append(thism4)

mass4 = np.array(mass4)

# plot mass above the density cuts as a function of time
ax.plot(times, mass1, color = c1, linewidth = 2)
ax.plot(times, mass2, color = c2, linewidth = 2)
#ax.plot(times, mass3, color = '0.6', linewidth = 2)
ax.plot(times, sinkmasses, color = c4, linewidth = 2)
horiz = alltlow + 0.05 * (allthi - alltlow)
vertbase = 1.5e4
ax.text(horiz, vertbase, r'n > $\mathdefault{10^3 cm^{-3}}$', transform=ax.transData,
    va = 'baseline', fontproperties = lfm, color=c1, snap = False) 
ax.text(horiz, vertbase / 2, r'n > $\mathdefault{10^4 cm^{-3}}$', transform=ax.transData,
    va = 'baseline', fontproperties = lfm, color=c2, snap = False)
ax.text(horiz, vertbase / 4, r'sinks', transform=ax.transData,
    va = 'baseline', fontproperties = lfm, color=c4, snap = False)             

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(alltlow, allthi)
ax.set_ylim(10.0, 3.e4)
set_ticks(ax, '0.6')

ax.set_xlabel('time / Myr', fontproperties = tfm, size = 15)
ax.set_ylabel('mass / '+r'M${_\odot}$', fontproperties = tfm, size = 15)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)

plt.savefig(outdir+'DenseGasStarMasses_logtime.pdf')
ax.set_xscale('linear')
plt.savefig(outdir+'DenseGasStarMasses_lintime.pdf')
plt.clf()




# plot q values
fig = plt.figure(figsize = (5,3.5))

ax2 = fig.add_axes([0.2, 0.2, 0.75, 0.75])

ax2.plot(times, np.array(mass1) / (np.array(sinkmasses) + 1.e-10), color = c1, linewidth = 2)
ax2.plot(times, np.array(mass2) / (np.array(sinkmasses) + 1.e-10), color = c2, linewidth = 2)
ax2.plot(times, np.array(mass4[:,0]) / (np.array(sinkmasses) + 1.e-10), color = c3, alpha = 0.7, linewidth = 2)
ax2.plot(times, np.array(mass4[:,1]) / (np.array(sinkmasses) + 1.e-10), color = c3, alpha = 0.7, linewidth = 2)
ax2.plot(times, np.array(mass4[:,2]) / (np.array(sinkmasses) + 1.e-10), color = c3, alpha = 0.7, linewidth = 2)
#ax2.plot(times, (np.array(mass1) + 0.7 * np.array(sinkmasses)) / (0.3 * np.array(sinkmasses)),
#    color = c1, linewidth = 0.75)
#ax2.plot(times, (np.array(mass2) + 0.7 * np.array(sinkmasses)) / (0.3 * np.array(sinkmasses)),
#    color = c2, linewidth = 0.75)
    
horiz = qtlow + 0.65 * (qthi - qtlow)    
ax2.text(horiz, 12.0, r'n > $\mathdefault{10^3 cm^{-3}}$', transform=ax2.transData,
    va = 'baseline', fontproperties = lfm, color=c1, snap = False) 
ax2.text(horiz, 10.0, r'n > $\mathdefault{10^4 cm^{-3}}$', transform=ax2.transData,
    va = 'baseline', fontproperties = lfm, color=c2, snap = False)
 
    
ax2.set_xscale('log')
ax2.set_xlim(qtlow, qthi)
ax2.set_ylim(0,15.0)
set_ticks(ax2, '0.6')

ax2.set_xlabel('time / Myr', fontproperties = tfm, size = 15)
ax2.set_ylabel('dense gas mass / sink mass', fontproperties = tfm, size = 15)

for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_fontproperties(tfm)

plt.savefig(outdir+'qRatio_logtime.pdf')
ax2.set_xscale('linear')
plt.savefig(outdir+'qRatio_lintime.pdf')



    
