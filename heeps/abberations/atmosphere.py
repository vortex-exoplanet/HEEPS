import numpy as np
import cv2
import proper
import astropy.units as u

def atmosphere(wfo, npupil, phase_screen):
    """ Compute the atmosphere phase screens as exp(iϕ(2π/λ)) and apply them to 
    the wavefront. The phase screen can be a single 2D array (nx, ny), or a cube 
    of N phase screens (N, nx, ny).
    """
    
    # pupil size must be odd (PROPER sets the center up-right next to the grid center)
    npupil = int(npupil + 1) if npupil % 2 == 0 else int(npupil)
    # reshape the phase screens to a cube [N × nx × ny] where N = 1 by default
    phase_screen = np.array(phase_screen, ndmin=3, dtype='float32')
    
    # get the number of screens
    ncube = phase_screen.shape[0]
    # get the numbers of pixels in x and y
    nx, ny = phase_screen.shape[1:3]
    # get the scaling factors from phase screens to simulation pupil size
    fx, fy = np.float64(npupil)/(nx, ny)
    # create a function for rescaling the phase screens
    rescale = lambda p: cv2.resize(p, (0,0), fx=fx, fy=fy, interpolation=cv2.INTER_LINEAR)
    # now, rescale the phase screens, and pad with zeros to match PROPER gridsize
    r = int((proper.prop_get_gridsize(wfo) - npupil)/2)
    phi = np.array([np.pad(rescale(phase_screen[x,:,:]), [(r+1,r),(r+1,r)], \
            mode='constant') for x in range(ncube)])
    
    # get the wavelength in microns
    lam = (proper.prop_get_wavelength(wfo)*u.m).to('um').value
    # get the wavenumber (spatial angular frequency)
    k = 2*np.pi/lam
    # get the complex atmosphere phase screens (TODO: allow cubes)
    atm_screen = np.exp(1j*k*phi[0])
    # multiply the atmosphere phase screens to the wavefront
    proper.prop_multiply(wfo, atm_screen)
    
    return wfo
