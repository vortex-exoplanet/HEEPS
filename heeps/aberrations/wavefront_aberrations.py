import heeps.util.img_processing as impro
import proper
import numpy as np
from astropy.io import fits
import astropy.units as u
from .static_ncpa import static_ncpa
import os.path

def wavefront_aberrations(wf, zern_inds=[2,3], zernike=None, npetals=6,
            petal_piston=None, atm_screen=None, ncpa_screen=None, **conf):
    
    # get useful parameters
    lam = conf['lam']
    gridsize = conf['gridsize']
    npupil = conf['npupil']
    
    """ tip-tilt = Zernike modes 2 and 3 """
    if zernike is not None:
        # translate the tip/tilt from lambda/D into RMS phase errors
        # RMS = x/2 = ((lam/D)*(D/2))/2 = lam/4
        zernike = np.array(zernike, ndmin=1)*lam/4
        zern_inds = np.array(zern_inds, ndmin=1)
        assert zern_inds.size == zernike.size, \
            "Zernike values and indices must be arrays/lists of same length."
        proper.prop_zernikes(wf, zern_inds, zernike)
    
    """ phase screen (atmosphere, petal piston) """
    if atm_screen is not None:
        atm_screen = np.nan_to_num(atm_screen)
        # resize and pad with zeros to match PROPER gridsize
        atm_screen = impro.resize_img(atm_screen, npupil)
        atm_screen = impro.pad_img(atm_screen, gridsize)
    else:
        atm_screen = 0
    if petal_piston is not None:
        # path to petal piston files
        filename = os.path.join(conf['input_dir'], '1024_pixelsize5mas_petal%s_243px.fits')
        # multiply all pistons by their respective petal, and sum them up
        piston_screen = np.sum(np.float32([petal_piston[x] \
                *fits.getdata(filename%(x+1)) for x in range(npetals)]), 0)
        # resize and pad with zeros to match PROPER gridsize
        piston_screen = impro.resize_img(piston_screen, npupil)
        piston_screen = impro.pad_img(piston_screen, gridsize)
    else:
        piston_screen = 0
    phase_screen = atm_screen + piston_screen
    if np.any(phase_screen) != 0:
        # wavenumber (spatial angular frequency) in rad / µm
        k = 2*np.pi/(lam*u.m).to('um').value
        # multiply the wavefront by the complex phase screen
        proper.prop_multiply(wf, np.exp(1j*k*phase_screen))
    
    """ static NCPA """
    if ncpa_screen is not None:
        ncpa_screen = np.nan_to_num(ncpa_screen)# in nm
        wf = static_ncpa(wf, ncpa_screen*1e-9, **conf)
    
    return wf
