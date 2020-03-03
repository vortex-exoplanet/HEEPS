import heeps.util.img_processing as impro
import proper
import numpy as np
from astropy.io import fits
import astropy.units as u
from heeps.aberrations.static_ncpa import static_ncpa
import os.path

def wavefront_aberrations(wf, zern_inds=[2,3], zernike=None, piston_screen=None, 
            atm_screen=None, ncpa_screen=None, resized=False, nans=True, **conf):
    
    # get useful parameters
    diam = conf['diam']
    gridsize = conf['gridsize']
    npupil = conf['npupil']
    
    """ tip-tilt = Zernike modes 2 and 3 """
    if zernike is not None:
        # convert the tip/tilt (mas) to RMS phase errors (m)
        # theta ~ tan(theta) = (2*RMS) / pupil radius => RMS = theta*D/4
        zernike = np.array(zernike, ndmin=1)*u.mas.to('rad')*diam/4
        zern_inds = np.array(zern_inds, ndmin=1)
        assert zern_inds.size == zernike.size, \
            "Zernike values and indices must be arrays/lists of same length."
        proper.prop_zernikes(wf, zern_inds, zernike)
    
    """ phase screen (atmosphere, petal piston) """
    if atm_screen is not None:
        if nans is True:
            atm_screen = np.nan_to_num(atm_screen)
        if resized is False:
            atm_screen = impro.resize_img(atm_screen, npupil)
    else:
        atm_screen = 0
    if piston_screen is not None:
        if nans is True:
            piston_screen = np.nan_to_num(piston_screen)
        if resized is False:
            piston_screen = impro.resize_img(piston_screen, npupil)
    else:
        piston_screen = 0
    if ncpa_screen is not None:
        if nans is True:
            ncpa_screen = np.nan_to_num(ncpa_screen)# in nm
        if resized is False:
            ncpa_screen = impro.resize_img(ncpa_screen, npupil)
    else:
        ncpa_screen = 0
    # sum resized phase screens
    phase_screen = atm_screen*1e-6 + piston_screen*1e-6 + ncpa_screen*1e-9
    
    if np.any(phase_screen) != 0:
        # pad with zeros to match PROPER gridsize
        phase_screen = impro.pad_img(phase_screen, gridsize)        
        # add the phase screen
        proper.prop_add_phase(wf, phase_screen)
    
    
    return wf
