from heeps.util.img_processing import resize_cube
from heeps.util.coord import mas2rms
import os.path
from astropy.io import fits
import numpy as np

def load_errors(nframes=20, nstep=1, npupil=285, add_phase=True, f_phase='',
        add_amp=False, f_amp='', add_point_err=False, f_point_err='', dit=0.3,
        add_apo_drift=False, apo_drift=0.02, apo_misalign=None, 
        add_ls_drift=False, ls_drift=0.02, ls_misalign=None,
        verbose=False, **conf):

    ''' Load wavefront errors.

    nframes (int):
        number of frames to crop the input data
    nstep (int):
        take 1 frame every nstep (cubesize = nframes/nstep)
    npupil (int):
        size of pupil
    add_phase (bool):
        true if adding phase screens
    f_phase (str):
        path to phase screens fits file
    add_amp (bool):
        true if adding amplitude screens
    f_amp (str):
        path to amplitude screens fits file
    add_point_err (bool):
        true if adding pointing errors
    f_point_err (str):
        path to pointing errors fits file
    dit (float):
        detector integration time in s
    add_apo_drift (bool):
        true if adding apodizer drift
    apo_drift (float):
        apodizer drift
    apo_misalign (list):
        constant apodizer misalignment (if no drift)
    add_ls_drift (bool):
        true if adding lyot stop drift
    ls_drift (float):
        lyot stop drift
    ls_misalign (list):
        constant lyot stop misalignment (if no drift)
    '''
    # size of cubes/arrays
    nscreens = int((nframes/nstep) + 0.5)

    # load phase screens
    if add_phase is True:
        assert(os.path.isfile(f_phase) and os.path.splitext(f_phase)[1] == '.fits'),\
            "'f_phase' must be a valid fits file."
        phase_screens = fits.getdata(f_phase)[:nframes][::nstep] # in meters
        phase_screens = np.nan_to_num(phase_screens)
        phase_screens = resize_cube(phase_screens, npupil)
        if verbose is True:
            print("Load phase screens from '%s'"%os.path.basename(f_phase))
            print('   nscreens=%s (nframes=%s, nstep=%s)'
                %(len(phase_screens), nframes, nstep))
    else:
        phase_screens = [None]*nscreens

    # load amp screens
    if add_amp is True:
        assert(os.path.isfile(f_amp) and os.path.splitext(f_amp)[1] == '.fits'),\
            "'f_amp' must be a valid fits file."
        amp_screens = np.array(fits.getdata(f_amp), ndmin=3)
        if len(amp_screens) > 1: # cube
            amp_screens = amp_screens[:nframes][::nstep]
        amp_screens = np.nan_to_num(amp_screens)
        amp_screens = resize_cube(amp_screens, npupil)
        if verbose is True:
            print("Load amp screens from '%s'"%os.path.basename(f_amp))
            print('   nscreens=%s'%(len(amp_screens)))
    else:
        amp_screens = np.array([None]*nscreens)

    # load pointing errors (in mas)
    if conf['disp'] is True:
        dispersion = os.path.join(conf['dir_temp'], 'dispersion_vals.fits')
        disp = np.array(fits.getdata(dispersion), ndmin= 2)
        if len(disp) > 1: # cube
            disp = disp[:nframes][::nstep]
        if verbose is True:
            print("Load dispersion values '%s'"%os.path.basename(dispersion))
            print('   nscreens=%s'%(len(disp)))
    else:
        disp = np.zeros((nscreens,2))#np.array([None]*nscreens)

    if add_point_err is True:
        tiptilts = np.array(fits.getdata(f_point_err), ndmin= 2)
        if len(tiptilts) > 1: # cube
            tiptilts = tiptilts[:nframes][::nstep]
        # convert mas to rms phase error
        tiptilts = mas2rms(tiptilts, conf['diam_ext'])
        if verbose is True:
            print("Load pointing errors from '%s'"%os.path.basename(f_point_err))
            print('   nscreens=%s'%(len(tiptilts)))
    else:
        tiptilts = np.zeros((nscreens,2))#np.array([None]*nscreens)

    # load apodizer drift
    if add_apo_drift is True and 'RAVC' in conf['mode']:
        apo_misaligns = np.array([[x,0,0,0,0,0] \
            for x in np.linspace(-apo_drift/2, apo_drift/2,
                                int(3600/dit))])[:nframes][::nstep]
        if verbose is True:
            print('Load apodizer drift=%s %% ptv'%apo_drift)
    else:
        apo_misaligns = np.array([apo_misalign]*nscreens)   # constant misalignment

    # load lyot stop drift
    if add_ls_drift is True and ('VC' in conf['mode'] or 'LC' in conf['mode']):
        ls_misaligns = np.array([[x,0,0,0,0,0] \
            for x in np.linspace(-ls_drift/2, ls_drift/2,
                                int(3600/dit))])[:nframes][::nstep]
        if verbose is True:
            print('Load Lyot stop drift=%s %% ptv'%ls_drift)
    else:
        ls_misaligns = np.array([ls_misalign]*nscreens)     # constant misalignment

    return phase_screens, amp_screens, tiptilts+disp, apo_misaligns, ls_misaligns