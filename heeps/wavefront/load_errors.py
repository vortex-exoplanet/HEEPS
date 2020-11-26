import os.path
from astropy.io import fits
import astropy.units as u
import numpy as np

def load_errors(nframes=20, nstep=1, npupil=285, add_phase=True, file_phase='', 
        add_point_err=False, file_point_err='', add_apo_drift=False, verbose=False, 
        ptv_drift=0.02, **conf):
    
    ''' Load wavefront errors.
    
    nframes: int
        number of frames to crop the input data
    nstep: int
        take 1 frame every nstep (cubesize = nframes/nstep)
    add_scao: bool
        true if adding SCAO residuals
    file_scao: str
        path to SCAO residuals fits file
    add_ncpa: bool
        true if adding NCPA (including petal piston)
    file_ncpa: str
        path to NCPA fits file
    add_pointing: bool
        true if adding pointing errors
    
    '''
    # size of cubes/arrays
    ncube = int((nframes/nstep) + 0.5)

    # load SCAO residuals
    if add_phase is True:
        assert(os.path.isfile(file_phase) and os.path.splitext(file_phase)[1] == '.fits'), \
            "'file_phase' must be a valid fits file."
        phase_screens = fits.getdata(file_phase)[:nframes][::nstep] # in meters
    else:
        phase_screens = [None]*ncube

    # load pointing errors
    if add_point_err is True:
        tiptilts = fits.getdata(file_point_err)[:nframes][::nstep]
    else:
        tiptilts = [None]*ncube

    # add apodizer drift
    if add_apo_drift is True:
        misaligns = np.array([[x,0,0,0,0,0] for x in np.linspace(-ptv_drift/2, ptv_drift/2, 12000)])[:nframes][::nstep]
    else:
        misaligns = [None]*ncube

    if verbose is True:
        print('Load wavefront error cubes of size %s (nframes=%s, nstep=%s)\n'\
            %(ncube, nframes, nstep))
    
    return phase_screens, tiptilts, misaligns