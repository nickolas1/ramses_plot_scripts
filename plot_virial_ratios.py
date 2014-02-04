from __future__ import division

import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import brewer2mpl
from os.path import expanduser
from matplotlib import cm
from astropy.io import ascii
         

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
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
    
def get_energies_and_virial_ratio(egyfile):
    """ reads in the reduced energy data and returns it, as well as the virial ratio.
        this is all in code units.
    """
    egy_data =  ascii.read(egyfile)
    epot = np.array(egy_data['epot'])
    ekin = np.array(egy_data['ekin'])
    eint = np.array(egy_data['eint'])    
    alphavir = 2*ekin/ np.abs(epot)
    return (epot, ekin, eint, alphavir)
    
def get_times(timefile, tff):
    """ reads in the reduced time data and returns it. this is in code units, 
        also returns a converted version to Myr and a version in tff.
        input tff is the freefall time in Myr.
    """
    time_data = ascii.read(timefile)
    times = np.array(time_data['time'])
    timesMyr = times * 14.909
    timesTff = timesMyr / tff
    return (times, timesMyr, timesTff)
    
    
def tff(density):
    """
    arguments: in density in Msun / pc^3
    returns: free-fall time in Myr
    """
    G = 6.67e-8
    pc = 3.086e18
    yr = 31.5576e6
    msun = 1.99e33
    denscgs = density * msun / pc**3
    tfreefall =  np.sqrt(3 * np.pi / (32 * G * denscgs))
    return tfreefall / yr / 1.e6 
    
    
def main():  
    fig = plt.figure(figsize = (5,3.5))
    ax = fig.add_axes([.18, .2, .75, .75])
    
    diffusedir = 'turbshock1024k4ghc7/'
    compactdir = 'turbshock1024k4ghcs7/'
    
    (epot, ekin, eint, alphavir) = get_energies_and_virial_ratio(diffusedir+'energies.log')
    fftime = tff(5.e4 / (4 / 3 * np.pi * 25**3))
    t100 = 4.0625 / fftime # time when 100 msun of stars have formed in Myr / tff
    t1000 = 5.12 / fftime
    (times, timesMyr, timesTff) = get_times(diffusedir+'mainsteptimes.log', fftime)
    
    (epot, ekin, eint, alphavir) = get_energies_and_virial_ratio(compactdir+'energies.log')
    fftime = tff(5.e4 / (4 / 3 * np.pi * 12.5**3))
    t100 = 1.1 / fftime
    t1000 = 1.53 / fftime
    (times, timesMyr, timesTff) = get_times(compactdir+'mainsteptimes.log', fftime)
    
    etot = epot + ekin + eint
    e0 = np.abs(etot[0])
    
    tplot = timesTff
    tmax = tplot[-1]
    ax.plot(tplot, np.abs(epot)/e0, color=c2, linewidth=1.75) #blue
    ax.plot(tplot, ekin/e0, color=c4, linewidth=1.75) #red
    # internal energy is boring, don't plot it.
    #ax.plot(tplot, eint/e0, color=c3, linewidth=1.5) #green
    ax.plot(tplot, alphavir, color=c1, linewidth=1.75)
    
    ax.set_ylim(.1, 10)
    ax.set_xlim(0, tmax)
    ax.set_yscale('log')
    ax.xaxis.grid(False,which='minor')
    #ax.yaxis.grid(False,which='minor')
    
    ax.set_xlabel(r'time / t$\mathdefault{_{ff}}$', fontproperties = tfm, size=15, color='0.15')
    ax.set_ylabel(r'$\alpha$, energy / E$\mathdefault{_{0}}$', fontproperties = tfm, size=15, color='0.15')
    
    # label the curves
    kinht = ekin[0] / e0
    ax.text(0.05*tmax, kinht, r'kinetic', transform=ax.transData,
         va = 'top', ha = 'left', fontproperties = lfm, color=c4, snap = False)
    potht = 1.05*(-epot[0] / e0)
    ax.text(0.05*tmax, potht, r'potential', transform=ax.transData,
         va = 'bottom', ha = 'left', fontproperties = lfm, color=c2, snap = False)
    alphaht = 1.3*alphavir[0]
    alphaht = alphavir[int(len(alphavir)/4)]
    ax.text(0.05*tmax, alphaht, r'virial ratio', transform=ax.transData,
         va = 'bottom', ha = 'left', fontproperties = lfm, color=c1, snap = False)
    
    # label the freefall time
    ax.text(0.7*tmax, 0.12, r't$\mathdefault{_{ff}}$ = %.1f Myr' % fftime,
        va = 'baseline', ha = 'left', fontproperties = lfm, color=c1, snap = False)
        
    # label 100 and 1000 Msun in sinks
    ax.plot([t100,t100],[.25, 7], color='0.5', zorder = 1)
    ax.plot([t1000,t1000],[.25, 7], color='0.5')
    
    set_ticks(ax, '0.15')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(tfm)
    
    plt.savefig('foo.png')
    
    sys.exit()
    epot = np.array(egy_data_diffuse['epot'])
    ekin = np.array(egy_data_diffuse['ekin'])
    eint = np.array(egy_data_diffuse['eint'])
    

    
     
        
if __name__ == '__main__':
    main()        