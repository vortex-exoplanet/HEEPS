import os.path
from astropy.io import fits
import astropy.units as u
import numpy as np

def load_errors(nframes=20, nstep=1, npupil=285, add_phase=True, add_amp=False, file_phase='', 
        file_amp='', add_point_err=False, file_point_err='', add_apo_drift=False, verbose=False, 
        ptv_drift=0.02, **conf):
    
    ''' Load wavefront errors.
    
    nframes: int
        number of frames to crop the input data
    nstep: int
        take 1 frame every nstep (cubesize = nframes/nstep)
    add_phase: bool
        true if adding phase screens
    file_phase: str
        path to phase screens fits file
    add_pointing: bool
        true if adding pointing errors
    
    '''
    # size of cubes/arrays
    ncube = int((nframes/nstep) + 0.5)

    # load phase screens
    if add_phase is True:
        assert(os.path.isfile(file_phase) and os.path.splitext(file_phase)[1] == '.fits'), \
            "'file_phase' must be a valid fits file."
        phase_screens = fits.getdata(file_phase)[:nframes][::nstep] # in meters
    else:
        phase_screens = [None]*ncube

    # load amp screens
    if add_amp is True:
        assert(os.path.isfile(file_amp) and os.path.splitext(file_amp)[1] == '.fits'), \
            "'file_amp' must be a valid fits file."
        amp_screens = fits.getdata(file_amp)
        if amp_screens.ndim == 3: # cube
            amp_screens = amp_screens[:nframes][::nstep]
    else:
        amp_screens = np.array([None]*ncube)
    
    # load pointing errors
    if add_point_err is True:
        tiptilts = fits.getdata(file_point_err)[:nframes][::nstep]
    else:
        tiptilts = np.array([None]*ncube)

    # add apodizer drift
    if add_apo_drift is True:
        misaligns = np.array([[x,0,0,0,0,0] \
            for x in np.linspace(-ptv_drift/2, ptv_drift/2, 12000)])[:nframes][::nstep]
    else:
        misaligns = np.array([None]*ncube)

    if verbose is True:
        print('Load wavefront error cubes of size %s (nframes=%s, nstep=%s)\n'\
            %(ncube, nframes, nstep))
    
    return phase_screens, amp_screens, tiptilts, misaligns