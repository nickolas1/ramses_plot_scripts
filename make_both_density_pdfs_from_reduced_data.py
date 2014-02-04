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
fontdir = homedir+'/Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=13)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=12)
    
fig = plt.figure(figsize = (10,3.5))

ax = fig.add_axes([0.1, 0.2, 0.35, 0.75])
ax2 = fig.add_axes([0.6, 0.2, 0.35, 0.75])

times = []
mass1 = []
mass2 = []
mass3 = []
sinkmasses = []


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

snapstart = int(sys.argv[1])
snapend = int(sys.argv[2])
snapiter = int(sys.argv[3])

for snap in range(snapstart, snapend, snapiter):
    datanamelo = 'reduced_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins_Low.dat'
    datanamemid = 'reduced_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins_Mid.dat'
    datanamehi = 'reduced_'+str(snap).zfill(5)+'/MassAndVolumeInDensityBins_High.dat'
    infoname = 'reduced_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'reduced_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    
    framesdir = outdir+'both_pdfs/'
    if not os.path.exists(framesdir):
        os.makedirs(framesdir)
    
    (time, unit_t) = get_time(infoname)
    
    # volume of the sphere the profile comes from
    (boxlen, unit_l) = get_boxsize(infoname)
    if boxlen > 7:
        sphererad = 27.0
    else:
        sphererad = 13.5
    spherevol = 4.0 * np.pi / 3.0 * (sphererad * 3.086e18)**3
    
    timeMyr = time * unit_t / 31557600.0 / 1.e6
    timeMyrRoundString = np.round(timeMyr)
    times.append(timeMyr)

    data = ascii.read(datanamelo)
    ax.plot(data['Lowdens'], data['CellVolume']/spherevol, color = colors[1], linewidth = 1.5)
    
    data = ascii.read(datanamemid)
    ax.plot(data['Middens'], data['CellVolume']/spherevol, color = colors[2], linewidth = 1.5)
    
    data = ascii.read(datanamehi)
    ax.plot(data['Highdens'], data['CellVolume']/spherevol, color = colors[3], linewidth = 1.5)
    
    # maximum likelihood fit to the first snap
    if snap == -1:
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
           

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1.e-24, 2.e-17)
    ax.set_ylim(1.e-9,0.1)
    set_ticks(ax, '0.6')

    ax.xaxis.grid(False,which='minor')
    ax.yaxis.grid(False,which='minor')

    ax.set_xlabel(r'density / g $\mathdefault{cm^{-3}}$', fontproperties = tfm, size = 15)
    ax.set_ylabel(r'volume weighted PDF', fontproperties = tfm, size = 15)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)

#plt.savefig(outdir+'VolumeDensityPDFs.pdf')

    mu = 2.33 # mean molecular weight
    mH = 1.6733e-24
    cdensmax = np.log10(10**2.0 / (mu * mH))  # convert these from g to n
    cdensmin = np.log10(10**-6.0 / (mu * mH))
    bins = 128
    binvals = np.arange(cdensmin, 1.000001*cdensmax, (cdensmax - cdensmin) / (bins))
    binmids = 0.5 * (np.roll(binvals, -1) + binvals)
    binmids = binmids[:len(binmids) - 1]
    
    files = [
        'surface_density_0.hdf5',
        'surface_density_1.hdf5',
        'surface_density_2.hdf5']
    files = ['surface_density_0.hdf5']
    colors = [c1,c2,c3,c4]
    
    f = h5py.File('reduced_'+str(snap).zfill(5)+'/surface_density_0.hdf5', 'r')
    sd = f['surface_density']
    totalhist = np.zeros(bins)
    # totalhistall keeps count of how many pixels are part of the cloud's area
    totalhistall = 0
    print snap
    for j in xrange(sd.shape[0]):
        coldensvals = sd[j]
        coldensvals -= np.log10(mu * mH)
        #if j == 600:
        #    print coldensvals
        #    print cdensmin,cdensmax
        hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
        # the maximum value of the histogram is the background ISM- don't count that 
        # in the calculation of the cloud area
        totalhistall += np.sum(hist) - np.max(hist)
        totalhist += hist  
    f.close() 
    #totalhistall = np.sum(totalhist)
    ax2.plot(10**binmids, totalhist/totalhistall, color = colors[0], linewidth = 1.5)
    print totalhistall
    
    # since we have three channels contributing to this PDF, need this factor of three
    # to correctly normalize them.
    totalhistall *= 3
    
    f = h5py.File('reduced_'+str(snap).zfill(5)+'/surface_density_low_0.hdf5', 'r')
    sd = f['surface_density']
    totalhist = np.zeros(bins)
    for j in xrange(sd.shape[0]):
        coldensvals = sd[j]
        coldensvals -= np.log10(mu * mH)
        #if j == 600:
        #    print coldensvals
        #    print cdensmin,cdensmax
        hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
        totalhist += hist  
    f.close() 
    ax2.plot(10**binmids, totalhist/totalhistall, color = colors[1], linewidth = 1.5)
    print np.sum(totalhist),totalhistall
   
    f = h5py.File('reduced_'+str(snap).zfill(5)+'/surface_density_mid_0.hdf5', 'r')
    sd = f['surface_density']
    totalhist = np.zeros(bins)
    for j in xrange(sd.shape[0]):
        coldensvals = sd[j]
        coldensvals -= np.log10(mu * mH)
        #if j == 600:
        #    print coldensvals
        #    print cdensmin,cdensmax
        hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
        totalhist += hist  
    f.close() 
    ax2.plot(10**binmids, totalhist/totalhistall, color = colors[2], linewidth = 1.5)
    print np.sum(totalhist),totalhistall
   
    f = h5py.File('reduced_'+str(snap).zfill(5)+'/surface_density_high_0.hdf5', 'r')
    sd = f['surface_density']
    totalhist = np.zeros(bins)
    for j in xrange(sd.shape[0]):
        coldensvals = sd[j]
        coldensvals -= np.log10(mu * mH)
        #if j == 600:
        #    print coldensvals
        #    print cdensmin,cdensmax
        hist, binedges = np.histogram(coldensvals, range = (cdensmin, cdensmax), bins = binvals)
        totalhist += hist  
    f.close() 
    ax2.plot(10**binmids, totalhist/totalhistall, color = colors[3], linewidth = 1.5)
    np.sum(totalhist),totalhistall
    print np.sum(totalhist),totalhistall
    
    ax2.set_xlabel(r'surface density / $\mathdefault{cm^{-2}}$', fontproperties = tfm, size = 15)
    #ax.set_ylabel('d', fontproperties = tfm, size = 15)
    ax2.set_ylabel(r'area weighted PDF', fontproperties = tfm, size = 15)
    
    ax2.set_xlim(0.99e19,1.01e24)
    ax2.set_ylim(0.99e-5,1.01e-1)
    ax2.set_yscale('log',nonposy='clip')
    ax2.set_xscale('log')
    set_ticks(ax2, '0.6')
    
    ax2.xaxis.grid(False,which='minor')
    ax2.yaxis.grid(False,which='minor')
    
    for label in ax2.get_xticklabels() + ax2.get_yticklabels():
        label.set_fontproperties(tfm)




    plt.savefig(framesdir+'BothDensityPDFs_'+str(snap).zfill(5)+'.png')
    plt.close()


    
