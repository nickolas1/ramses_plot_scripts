from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
from astropy.io import ascii
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'/Documents/astronomy/macfontsforpython/'
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

# volume of the sphere the profile comes from
sphererad = 27.0
spherevol = 4.0 * np.pi / 3.0 * (sphererad * 3.086e18)**3

print os.path.basename(os.path.normpath(os.getcwd()))
if os.path.basename(os.path.normpath(os.getcwd())) == 'turbshock512k4gcl':
    snaps = [15, 82, 95, 109]
    myrstrings = ['1', '6', '7', '8']
if os.path.basename(os.path.normpath(os.getcwd())) == 'turbshock512k4gcsl':
    snaps = [4, 28, 41, 55]
    myrstrings = ['0.1', '1.0', '1.5', '2.0']
colors = [c1, c2, c3, c4]

#snaps=[1]
#myrstrings=['0.0']

for i in range(len(snaps)):
    snap = snaps[i]
    dataname = 'output_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins.dat'
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    
    if os.path.isfile(dataname):
        (time, unit_t) = get_time(infoname)
        timeMyr = time * unit_t / 31557600.0 / 1.e6
        timeMyrRoundString = np.round(timeMyr)
        times.append(timeMyr)
    
        data = ascii.read(dataname)
        
        ax.plot(data['Density'], data['CellVolume']/spherevol, color = colors[i], linewidth = 1.5)
        
        # maximum likelihood fit to the first snap
        if i == 1:
            dens = np.array(data['Density'])
            vals = np.array(data['CellVolume'])
            # find the maximum value above 1.e-24 (want to avoid the ISM peak)
            cloudsel = (dens > 1.e-24)
            denssub = dens[cloudsel]
            valssub = vals[cloudsel]
            densmax = denssub[valssub.argmax()]
            logdensmax = np.log10(densmax)

            # select 3 decades surrounding the density maximum to fit to
            logdens = np.log10(dens)
            fitsel = logdens[np.abs(logdens - logdensmax) <= 1.5]
            print fitsel
            # get the maximum likelihood lognormal to those 3 decades
            muML = fitsel.sum() / len(fitsel)
            print len(fitsel)
            print muML

            sys.exit()
            # fit to three orders of magnitude surrounding peak density
            maxpt = data['Density'].argmax()
            maxdens = np.max(data['Density'])
            fitvals = np.max
        
horiz = 1.2e-24
vertbase = 1.e-5
ax.text(horiz, vertbase, myrstrings[0]+' Myr', transform=ax.transData,
    va = 'baseline', fontproperties = lfm, color=c1, snap = False) 
if len(myrstrings) > 1:
	ax.text(horiz, vertbase / 5, myrstrings[1]+' Myr', transform=ax.transData,
    	    va = 'baseline', fontproperties = lfm, color=c2, snap = False)
	ax.text(horiz, vertbase / 5**2, myrstrings[2]+' Myr', transform=ax.transData,
    	    va = 'baseline', fontproperties = lfm, color=c3, snap = False)      
	ax.text(horiz, vertbase / 5**3, myrstrings[3]+' Myr', transform=ax.transData,
    	    va = 'baseline', fontproperties = lfm, color=c4, snap = False)          

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1.e-25, 2.e-18)
ax.set_ylim(1.e-8,0.1)
set_ticks(ax, '0.6')

ax.xaxis.grid(False,which='minor')
ax.yaxis.grid(False,which='minor')

ax.set_xlabel(r'density / g $\mathdefault{cm^{-3}}$', fontproperties = tfm, size = 15)
#ax.set_ylabel('d', fontproperties = tfm, size = 15)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontproperties(tfm)

plt.savefig(outdir+'VolumeDensityPDFs.pdf')




    
