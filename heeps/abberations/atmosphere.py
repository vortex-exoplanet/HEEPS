import numpy as np
import proper
import astropy.units as u
from skimage.transform import resize

def atmosphere(wfo, npupil, phase_screen):
    """ Compute the atmosphere phase screens as exp(iϕ(2π/λ)) and apply them to 
    the wavefront. The phase screen can be a single 2D array (nx, ny), or a cube 
    of N phase screens (N, nx, ny).
    """
    
    # check pupil size is odd (PROPER sets the center up-right next to the grid center)
    npupil = int(npupil + 1) if npupil % 2 == 0 else int(npupil)
    # check phase screen shape is 3D [N × nx × ny] where N = 1 by default
    phase_screen = np.array(phase_screen, ndmin=3, dtype='float32')
    
    # TODO: allow multiple wavelengths if ncube = (1 or size(lam))
    
    # get the wavelength in microns
    lam = (proper.prop_get_wavelength(wfo)*u.m).to('um').value
    # get the wavenumber (spatial angular frequency)
    k = 2*np.pi/lam
    # get the number of screens
    ncube = phase_screen.shape[0]
    # create a function for rescaling the phase screens
    rescale = lambda p: resize(p, (npupil, npupil), order=1, preserve_range=True)
    
    # for each phase screen:
    # - rescale to pupil size, 
    # - pad with zeros to match PROPER gridsize
    # - multiply the wavefront by the complex atmosphere phase screen
    r = int((proper.prop_get_gridsize(wfo) - npupil)/2)
    phi = np.array([np.pad(rescale(phase_screen[x,:,:]), [(r+1,r),(r+1,r)], \
            mode='constant') for x in range(ncube)])
    proper.prop_multiply(wfo, np.exp(1j*k*phi[0]))
    
    # TODO: allow cubes instead of using phi[0] (PROPER doesn't allow multiple wavefronts)
    
    return wfo
