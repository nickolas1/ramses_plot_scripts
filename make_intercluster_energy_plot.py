from __future__ import division

from yt.mods import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gc
import sys
from os.path import expanduser

# import ramses helper functions and get figure directory
homedir = '/home/moon/moeckel/'
homedir = expanduser('~')+'/'
print homedir
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *
#mpl.rc_file(homedir+'pythonhelpers/ramses/matplotlibrc')
outdir = get_output_path(homedir)

# set some fonts
fontdir = homedir+'Documents/astronomy/macfontsforpython/'
tfm = fm.FontProperties( # tick font main
    fname=fontdir+'Gotham-Book.ttf', size=13)    
lfm = fm.FontProperties( # label font main
    fname=fontdir+'Gotham-BookItalic.ttf', size=11)  

# choose dark or light background
background = 'light'

if background == 'dark':
    imshowmap = 'gray'
    textcolor = c3l
    textalpha = 0.5
if background == 'light':
    imshowmap = 'gray_r'
    textcolor = textlightbg
    textalpha = 0.7
    
def get_total_energies(sinks):
    grav = 6.67e-8
    nsink = len(sinks)
    pot = 0.0
    kin = 0.5 * sinks[nsink - 1, 1] * np.sum(sinks[nsink - 1, 5:8]**2)
    for i in xrange(nsink - 1):
        imass = sinks[i, 1]
        ipos = sinks[i, 2:5]
        ivel = sinks[i, 5:8]
        kin += 0.5 * imass * np.sum(ivel**2)      
        for j in xrange(i+1, nsink):
            jmass = sinks[j, 1]
            jpos = sinks[j, 2:5]
            dr = np.sqrt(np.sum((ipos - jpos)**2))
            pot -= grav * imass * jmass / dr
    return(kin, pot)

def get_individual_energies(sinks):
    grav = 6.67e-8
    nsink = len(sinks)
    pots = 0.0
    kins = 0.5 * sinks[nsink - 1, 1] * np.sum(sinks[nsink - 1, 5:8]**2)
    for i in xrange(nsink - 1):
        imass = sinks[i, 1]
        ipos = sinks[i, 2:5]
        ivel = sinks[i, 5:8]
        kin += 0.5 * imass * np.sum(ivel**2)      
        for j in xrange(i+1, nsink):
            jmass = sinks[j, 1]
            jpos = sinks[j, 2:5]
            dr = np.sqrt(np.sum((ipos - jpos)**2))
            pot -= grav * imass * jmass / dr
    return(kin, pot)
    
fig = plt.figure(figsize = (9, 3))
ax1 = fig.add_axes([.05, .05, .27, .9])
ax2 = fig.add_axes([.37, .05, .27, .9])
ax3 = fig.add_axes([.69, .05, .27, .9])

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    framename = outdir+'framesD/frame_'+str(snap).zfill(4)+'.png'
    
    # see if we have any sink particles
    try:
        with open(sinkname): 
            sinks = get_sinks(sinkname)
            if len(sinks) > 0:                
                # move everything to center of momentum
                rcm = np.sum(sinks[:,1,np.newaxis] * sinks[:,2:5], axis=0) / np.sum(sinks[:,1])
                vcm = np.sum(sinks[:,1,np.newaxis] * sinks[:,5:8], axis=0) / np.sum(sinks[:,1])         
                sinks[:,2:5] -= rcm
                sinks[:,5:8] -= vcm
                
                # convert everything to cgs
                # masses are in Msun, distances in 10pc, vels in kms
                sinks[:,1] *= 1.99e33
                sinks[:,2:5] *= 10 * 3.086e18
                sinks[:,5:8] *= 1.e5
                
                (totalkin, totalpot) = get_total_energies(sinks)
                print totalkin, totalpot, totalkin / np.abs(totalpot)
                ax1.scatter(sinks[:,2],sinks[:,3],marker='.',s=5,facecolor=csink,edgecolor=csink)
                ax2.scatter(sinks[:,2],sinks[:,4],marker='.',s=5,facecolor=csink,edgecolor=csink)
                ax3.scatter(sinks[:,3],sinks[:,4],marker='.',s=5,facecolor=csink,edgecolor=csink)
                
                # add one Myr tracks
                sinks2 = sinks.copy()
                Myrs = 1
                sinks2[:,2:5] += sinks[:,5:8] * Myrs * 1.e6 * 31557600.0
                for i in xrange(len(sinks)):
                    ax1.arrow(sinks[i,2],sinks[i,3],sinks2[i,2]-sinks[i,2],sinks2[i,3]-sinks[i,3],lw=0.4,head_width=3.e17,length_includes_head=True)
                    ax2.arrow(sinks[i,2],sinks[i,4],sinks2[i,2]-sinks[i,2],sinks2[i,4]-sinks[i,4],lw=0.4,head_width=3.e17,length_includes_head=True)
                    ax3.arrow(sinks[i,3],sinks[i,4],sinks2[i,3]-sinks[i,3],sinks2[i,4]-sinks[i,4],lw=0.4,head_width=3.e17,length_includes_head=True)
    except IOError:
        pass      
    
    ax1.set_xlim(-3.e19,3.e19)
    ax1.set_ylim(-3.e19,3.e19)    
    ax2.set_xlim(-3.e19,3.e19)
    ax2.set_ylim(-3.e19,3.e19)    
    ax3.set_xlim(-3.e19,3.e19)
    ax3.set_ylim(-3.e19,3.e19)    
    plt.savefig(outdir+'sinksfoo.png', dpi = 200)
        

    
