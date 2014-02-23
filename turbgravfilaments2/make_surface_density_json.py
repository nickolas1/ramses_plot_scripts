from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.cm as cm
import gc
import sys
import h5py
from astropy import constants as const
from astropy import units as u
from os.path import expanduser
from matplotlib import rcParams
import scipy.ndimage as ndimage


"""
plot column density from different density ranges. 

the hdf5 files contain log surface density in g cm^-2
files with 'C18O': surface density of gas with number density 10^3 < n < 10^4.5
files with 'N2Hplus': surface density of gas with number density 10^4.5 < n
files that just say 'surface_density': all the gas

this script plots column density in cm^-2, converting from g cm^-2 using mu=2.33
"""


# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *

#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font
    fname=fontdir+'Gotham-Book.ttf', size=7)
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=6) 

rcParams['xtick.direction'] = 'out'

outdir = get_output_path(homedir)
outdir = './'

snapstr = str(int(sys.argv[1])).zfill(5)
infoname = 'reduced_'+snapstr+'/info_'+snapstr+'.txt'
(boxlen, unit_l) = get_boxsize(infoname)


"""
   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   toggle these to plot a colorbar and length scale bar
   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
"""
PlotColorBars = False
PlotScaleBar = False
PlotSinks = False

outputres = 256

mu = 2.33

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    # first do projections
    imshowmap = 'nickmapVD2'
    #imshowmap = 'bone_r'
    
    cdmin = 10**-3.3 
    cdmax = 10**-1.5
    cdmin = 0
    cdmax = 4
    
    cdmin = 0
    cdmax = 1.*10**22
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'
    sinkname = fileprefix+'sink_'+str(snap).zfill(5)+'.csv'
    
    
    # downsampled surface density
    dsSD = np.zeros([outputres, outputres])
    # downsampled C18O
    dsC18O = np.zeros([outputres, outputres])
    # downsampled N2Hplus
    dsN2Hplus = np.zeros([outputres, outputres])
    
    for i in xrange(2,3):
        fileSD = fileprefix+'surface_density'+str(i)+'.hdf5'
        fileCO = fileprefix+'surface_density_C18O'+str(i)+'.hdf5'
        fileN2Hplus = fileprefix+'surface_density_N2Hplus'+str(i)+'.hdf5'
        if os.path.exists(fileCO):
            print snap,fileCO
            fSD = h5py.File(fileSD, 'r')
            sd = fSD['surface_density']
            # convert to linear units, divide by mu * mH
            sd = 10**np.array(sd) / (mu * const.m_p.cgs.value)
            sd = sd.transpose()
            # there are a lot of very small values- cull them before getting some
            # info on the range of interesting values  
            sdnonzero = sd[sd > 10**4]
            print 'mean non-zero column density: ',np.mean(sdnonzero),'cm^-2'
            print 'median non-zero column density: ',np.median(sdnonzero),'cm^-2'
            print 'max column density: ',np.max(sd),'cm^-2'
            
            fC18O = h5py.File(fileCO, 'r')
            sdC18O = fC18O['surface_density_C18O']
            # convert to linear units, divide by mu * mH
            sdC18O = 10**np.array(sdC18O) / (mu * const.m_p.cgs.value)
            sdC18O = sdC18O.transpose()
            # there are a lot of very small values- cull them before getting some
            # info on the range of interesting values  
            sdnonzero = sdC18O[sdC18O > 10**4]
            print 'mean non-zero C18O column density: ',np.mean(sdnonzero),'cm^-2'
            print 'median non-zero C18O column density: ',np.median(sdnonzero),'cm^-2'
            print 'max C18O column density: ',np.max(sdC18O),'cm^-2'
            
            fN2Hplus = h5py.File(fileN2Hplus, 'r')
            sdN2Hplus = fN2Hplus['surface_density_N2Hplus']
            # convert to linear units, divide by mu * mH
            sdN2Hplus = 10**np.array(sdN2Hplus) / (mu * const.m_p.cgs.value)
            sdN2Hplus = sdN2Hplus.transpose()
            # there are a lot of very small values- cull them before getting some
            # info on the range of interesting values  
            sdnonzero = sdN2Hplus[sdN2Hplus > 10**4]
            print ''
            print 'mean non-zero N2Hplus column density: ',np.mean(sdnonzero),'cm^-2'
            print 'median non-zero N2Hplus column density: ',np.median(sdnonzero),'cm^-2'
            print 'max N2Hplus column density: ',np.max(sdN2Hplus),'cm^-2'
            
            glom = sd.shape[0] / outputres
            for ii in xrange(dsSD.shape[0]):
                iindexlo = ii * glom
                iindexhi = (ii + 1) * glom
                for jj in xrange(dsSD.shape[1]):
                    jindexlo = jj * glom
                    jindexhi = (jj +1) * glom
                    dsSD[ii, jj] = sd[iindexlo:iindexhi, jindexlo:jindexhi].sum() / glom**2
                    dsC18O[ii, jj] = sdC18O[iindexlo:iindexhi, jindexlo:jindexhi].sum() / glom**2
                    dsN2Hplus[ii, jj] = sdN2Hplus[iindexlo:iindexhi, jindexlo:jindexhi].sum() / glom**2
            
            # convert 0 to 10**-99
            dsSD[dsSD < 1.e-99] = 1.e-99
            dsC18O[dsC18O < 1.e-99] = 1.e-99
            dsN2Hplus[dsN2Hplus < 1.e-99] = 1.e-99
            
            dsSD = np.round(np.log10(dsSD),2)
            dsC18O = np.round(np.log10(dsC18O),2)
            dsN2Hplus = np.round(np.log10(dsN2Hplus),2)
            
            # convert 10**-99 to 0
            dsSD[dsSD < 1] = 0
            dsC18O[dsC18O < 1] = 0
            dsN2Hplus[dsN2Hplus < 1] = 0
            
            f = open('columndensities_'+str(snap).zfill(5)+'.json', 'w')
            f.write('[\n')
            for ii in xrange(dsSD.shape[0]):
                f.write('{"x":'+str(ii)+',\n"SD":['+str(dsSD[ii, 0]))
                for jj in xrange(1, dsSD.shape[1]):
                    f.write(','+str(dsSD[ii, jj]))
                f.write('],\n')
                f.write('"C18O":['+str(dsC18O[ii, 0]))
                for jj in xrange(1, dsC18O.shape[1]):
                    f.write(','+str(dsC18O[ii, jj]))
                f.write('],\n')
                f.write('"N2Hplus":['+str(dsN2Hplus[ii, 0]))
                for jj in xrange(1, dsN2Hplus.shape[1]):
                    f.write(','+str(dsN2Hplus[ii, jj]))
                f.write(']\n')
                f.write('}')
                if ii < dsSD.shape[0] - 1:
                    f.write(',\n')
            f.write('\n]')
            f.close()
            