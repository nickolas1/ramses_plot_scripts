from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
from astropy.io import ascii

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
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)
    
fig = plt.figure(figsize = (5,3.5))

ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
#ax = fig.add_axes([0., 0., 1., 1.])

times = []
mass1 = []
mass2 = []
mass3 = []
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
    dataname = 'output_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins.dat'
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
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
        timeMyr = time * unit_t / 31557600.0 / 1.e6
        times.append(timeMyr)
    
        data = ascii.read(dataname)
    
        # first get M(>rho) / tff(rho)  
        cumulativemass = np.cumsum(data['CellMassMsun'])
        massabove = cumulativemass.max() - cumulativemass
        logdens = np.log10(data['Density'])
        G = 6.67e-8
        yr = 31557600.0
        msun = 1.99e33
        tffcgs =  np.sqrt(3 * np.pi / (32 * G * data['Density']))
        tffMyr = tffcgs / 1.e6 / yr 
        ndense = data['Density'] / (mu * mH)
        print massabove
        
        # add the mass in sinks to these massabove values
        sinks = get_sinks(sinkname)
        Msinks = np.sum(sinks[:,1])
        massabove[:] += Msinks
        rhoTffRatio = massabove / tffMyr
        
        # next approximate the instantaneous star formation rate
        infolast = 'output_'+str(snap-1).zfill(5)+'/info_'+str(snap-1).zfill(5)+'.txt'
        infonext = 'output_'+str(snap+1).zfill(5)+'/info_'+str(snap+1).zfill(5)+'.txt'
        sinklast = 'output_'+str(snap-1).zfill(5)+'/sink_'+str(snap-1).zfill(5)+'.out'
        sinknext = 'output_'+str(snap+1).zfill(5)+'/sink_'+str(snap+1).zfill(5)+'.out'
        sinkslast = get_sinks(sinklast)
        sinksnext = get_sinks(sinknext)
        Msinkslast = np.sum(sinkslast[:,1])
        Msinksnext = np.sum(sinksnext[:,1])
        (time, unit_t) = get_time(infolast)
        Tlast = time * unit_t / yr / 1.e6
        (time, unit_t) = get_time(infonext)
        Tnext = time * unit_t / yr / 1.e6
        sfeNow = (Msinksnext - Msinkslast) / (Tnext - Tlast)
        print Msinksnext, Msinkslast, Tnext, Tlast, sfeNow
        
        
ax.plot(ndense, sfeNow / rhoTffRatio, color = c2, linewidth = 2)

(time, unit_t) = get_time(infoname)
timeMyr = time * unit_t / yr / 1.e6
horiz = 1.e4
vert = 0.8
ax.text(horiz, vert, r'%.1f' %timeMyr, transform = ax.transData, 
    ha = 'right',va = 'baseline', fontproperties = lfm, color = c1, snap = False)
ax.text(1.1*horiz, vert, r'Myr', transform = ax.transData,
    ha = 'left', va = 'baseline', fontproperties = lfm, color = c1, snap = False)     
ax.set_ylim(.001, 1.0)
ax.set_xscale('log')
ax.set_yscale('log')        
    
    
set_ticks(ax, '0.6')

ax.set_xlabel('n / cm'+r'$\mathdefault{{-3}}$', fontproperties = tfm, size = 15)
ax.set_ylabel(r'$\mathdefault{SFR_{ff}}$', fontproperties = tfm, size = 15)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)

plt.savefig(outdir+'SFRff.pdf')


sys.exit()
 



    
