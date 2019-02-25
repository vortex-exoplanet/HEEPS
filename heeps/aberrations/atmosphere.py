import heeps.util.img_processing as impro
import proper
import numpy as np
import astropy.units as u

def atmosphere(wf, atm_screen, **conf):
    """ Compute the atmosphere phase screens as exp(iϕ(2π/λ)) and apply them to 
    the wavefront. 
    """
    
    # get useful parameters
    lam = conf['lam']
    gridsize = conf['gridsize']
    npupil = conf['npupil']
    ncube = atm_screen.shape[0]
    
    # resize and pad with zeros to match PROPER gridsize
    atm_screen = impro.resize_img(atm_screen, npupil)
    atm_screen = impro.pad_img(atm_screen, gridsize)
    
    # wavenumber (spatial angular frequency) in rad / µm
    k = 2*np.pi/(lam*u.m).to('um').value
    # multiply the wavefront by the complex phase screen
    proper.prop_multiply(wf, np.exp(1j*k*atm_screen))
    
    return wf
