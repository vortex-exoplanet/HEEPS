import os.path
from astropy.io import fits
import astropy.units as u

def load_errors(nframes=20, nstep=1, add_scao=False, file_scao='', 
        add_ncpa=False, file_ncpa='', add_petal_piston=False, file_petal_piston='',
        add_pointing_err=False, add_apo_drift=False, verbose=False, **conf):
    
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
    if add_scao is True:
        assert(os.path.isfile(file_scao) and os.path.splitext(file_scao)[1] == '.fits'), \
            "'file_scao' must be a valid fits file."
        scao_screens = fits.getdata(file_scao)[:nframes][::nstep]
    else:
        scao_screens = [None]*ncube

    # load NCPA
    if add_ncpa is True:
        assert(os.path.isfile(file_ncpa) and os.path.splitext(file_ncpa)[1] == '.fits'), \
            "'file_ncpa' must be a valid fits file."
        ncpa_screens = fits.getdata(file_ncpa)[:nframes][::nstep]
    else:
        ncpa_screens = [None]*ncube

    # load petal piston
    if add_petal_piston is True:
        assert(os.path.isfile(file_petal_piston) and os.path.splitext(file_petal_piston)[1] == '.fits'), \
            "'file_petal_piston' must be a valid fits file."
        petal_piston_screens = fits.getdata(file_petal_piston)[:nframes][::nstep]
    else:
        petal_piston_screens = [None]*ncube
    
    # TODO: 
    # combine phase screens
    # add pointing errors
    # add apodizer drift
    phase_screens = scao_screens
    pointing_errs = [None]*ncube
    apo_drifts = [None]*ncube

    if verbose is True:
        print('Load wavefront error cubes of size %s (nframes=%s, nstep=%s)\n'\
            %(ncube, nframes, nstep))
    
    return phase_screens, pointing_errs, apo_drifts