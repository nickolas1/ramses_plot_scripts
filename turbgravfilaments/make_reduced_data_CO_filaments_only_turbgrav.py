from __future__ import division

from yt.mods import *
from yt.config import ytcfg
import gc
import sys
import h5py
import shutil
from os.path import expanduser
from mpi4py import MPI

# import ramses helper functions and get figure directory
homedir = expanduser('~')+'/'

# import ramses helper functions and get figure directory
sys.path.append(homedir+'pythonhelpers/ramses/')
from ramses_helpers import *

def _CO(field, data):
    mu = 2.33 # mean molecular weight
    mH = 1.6733e-24
    lolim = .001 * mu * mH # not interested in anything below 10^3 / cm^3
    hilim = 31622.0 * mu * mH # not interested in anything above 10^4.5 / com^3
    newfield = data['Density']
    antiselection = (data['Density'] < lolim) | (data['Density'] >= hilim)
    newfield[antiselection] = 1.e-99
    
    # cut out the filaments we're interested in
    leftpoint = np.array([0.795, 0.34])
    rightpoint = np.array([0.805, 0.05]) 
    width = 0.045
    vector = rightpoint - leftpoint
    midpoint = leftpoint + 0.5 * vector
    # translate to midpoint
    transx = data['x'] - midpoint[0]
    transy = data['y'] - midpoint[1]
    length = np.linalg.norm(vector)
    orthovec = (-vector[1], vector[0])
    orthovec /= np.linalg.norm(orthovec)
    vector /= np.linalg.norm(vector)
    
    # rotate around midpoint. orthovec is already a unit vector now.
    beta = np.arccos(orthovec[1])
    rotx = transx * np.cos(beta) - transy * np.sin(beta)
    roty = transx * np.sin(beta) + transy * np.cos(beta)
    
    # cut based on width and length of box
    antiselection2 = (np.abs(rotx) > 0.5*length) | (np.abs(roty) > 0.5 * width)
    newfield[antiselection2] = 1.e-99
    
    return newfield

for snap in range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
    infoname = 'output_'+str(snap).zfill(5)+'/info_'+str(snap).zfill(5)+'.txt'
    sinkname = 'output_'+str(snap).zfill(5)+'/sink_'+str(snap).zfill(5)+'.out'
    
    fileprefix = 'reduced_'+str(snap).zfill(5)+'/'

    if ytcfg.getint('yt', '__topcomm_parallel_rank') == 0:
        if not os.path.exists(fileprefix):
            os.makedirs(fileprefix)       
            # copy the infofile and sinkfile to the reduced directory 
            shutil.copy(infoname, fileprefix)
            if os.path.exists(sinkname):
                shutil.copy(sinkname, fileprefix)
    
    (lmin, lmax) = get_level_min_max(infoname)
    (boxlen, unit_l) = get_boxsize(infoname)

    ds = load(infoname, fields=['Density','x-velocity','y-velocity','z-velocity','Pressure'])
    
    add_field('CO', function=_CO)
    
    # center on original center of cloud
    cntr = [0.5, 0.5, 0.5]

    wd = 1.0
    # res should be base resolution times 2**levels of refinement * wd
    lmaxplot = min(11, lmax)  
    resx = int(wd * 2**lmaxplot)
    expandfac = 1
    res = (expandfac * resx,expandfac * resx)
    ht = wd * res[0] / res[1]
    width = (wd, 'unitary')
    height = (ht, 'unitary')

    # cut out the filaments we're interested in
    leftpoint = np.array([0.795, 0.34])
    rightpoint = np.array([0.805, 0.05]) 
    width = 0.045
    vector = rightpoint - leftpoint
    midpoint = leftpoint + 0.5 * vector
    orthovec = (-vector[1], vector[0], 0)
    orthovec /= np.linalg.norm(orthovec)
    print orthovec
    cntr = [midpoint[0], midpoint[1], 0.5]
    width = [length, length, 1.0]
    #get a projection orthogonal to the filament
    proj = off_axis_projection(ds,
        cntr, #center
        orthovec, #normal vector
        width, #width
        256, #resolution
        'CO',
        north_vector = [0, 0, -1] #make z the vertical direction on the plot
        )
    print proj.shape
    """proj = off_axis_projection(ds,
        [0.5, 0.5, 0.5], #center
        orthovec, #normal vector
        1.0, #width
        res, #resolution
        'CO',
        north_vector = [0, 0, 1] #make z the vertical direction on the plot
        )"""
    proj.write_hdf5(fileprefix+'surface_density_CO_fil0.hdf5')

    del(proj)   
    del(ds)
    gc.collect()
    
