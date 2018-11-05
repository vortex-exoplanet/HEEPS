import numpy as np
import proper
import astropy.units as u
from skimage.transform import resize

def atmosphere(wfo, phase_screen, **conf):
    """ Compute the atmosphere phase screens as exp(iϕ(2π/λ)) and apply them to 
    the wavefront. The phase screen can be a single 2D array (nx, ny), or a cube 
    of N phase screens (N, nx, ny).
    """
    
    # TODO: allow multiple wavelengths if ncube = (1 or size(lam))
    # TODO: allow cubes instead of using phi[0] (PROPER doesn't allow multiple wavefronts)
    
    # check phase screen shape is 3D [N × nx × ny] where N = 1 by default
    phase_screen = np.array(phase_screen, ndmin=3, dtype='float32')
    
    # load useful parameters
    npupil = conf['npupil']
    lam = (conf['WAVELENGTH']*u.m).to('um').value
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
    
    return wfo
