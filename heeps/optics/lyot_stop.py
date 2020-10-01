import heeps.util.img_processing as impro
import numpy as np
import proper
from astropy.io import fits
import os.path

def lyot_stop(wf, mode='RAVC', ravc_r=0.6, ls_dRext=0.03, ls_dRint=0.05, 
        ls_dRspi=0.04, spi_width=0.5, spi_angles=[0,60,120], diam_ext=37, 
        diam_int=11, ls_misalign=None, file_app_phase='', file_app_amp='', 
        ngrid=1024, npupil=285, margin=50, get_amp=False, 
        get_phase=False, verbose=False, **conf):

    """Add a Lyot stop, or an APP."""
    
    # case 1: Lyot stop
    if mode in ['CVC', 'RAVC']:
        # LS parameters
        r_obstr = ravc_r if mode in ['RAVC'] else diam_int/diam_ext
        ls_int = r_obstr + ls_dRint
        ls_ext = 1 - ls_dRext
        ls_spi = spi_width/diam_ext + ls_dRspi
        # LS misalignments
        ls_misalign = [0,0,0,0,0,0] if ls_misalign is None else list(ls_misalign)
        dx_amp, dy_amp, dz_amp = ls_misalign[0:3]
        dx_phase, dy_phase, dz_phase = ls_misalign[3:6]
        # create Lyot stop
        proper.prop_circular_aperture(wf, ls_ext, dx_amp, dy_amp, NORM=True)
        if diam_int > 0:
            proper.prop_circular_obscuration(wf, ls_int, dx_amp, dy_amp, NORM=True)
        if spi_width > 0:
            for angle in spi_angles:
                proper.prop_rectangular_obscuration(wf, ls_spi, 2, \
                        dx_amp, dy_amp, ROTATION=angle, NORM=True)
        if verbose is True:
            print('Create Lyot stop')
            print('   ls_int=%3.4f, ls_ext=%3.4f, ls_spi=%3.4f'\
                %(ls_int, ls_ext, ls_spi))
            print('')

    # case 2: APP
    elif mode in ['APP']:
        if verbose is True:
            print('Load APP from files\n')
        # get amplitude and phase data
        APP_amp = fits.getdata(file_app_amp) if os.path.isfile(file_app_amp) \
                else np.ones((npupil, npupil))
        APP_phase = fits.getdata(file_app_phase) if os.path.isfile(file_app_phase) \
                else np.zeros((npupil, npupil))
        # resize to npupil
        APP_amp = impro.resize_img(APP_amp, npupil)
        APP_phase = impro.resize_img(APP_phase, npupil)
        # pad with zeros to match PROPER ngrid
        APP_amp = impro.pad_img(APP_amp, ngrid, 1)
        APP_phase = impro.pad_img(APP_phase, ngrid, 0)
        # multiply the loaded APP
        proper.prop_multiply(wf, APP_amp*np.exp(1j*APP_phase))
    
    # get the LS amplitude and phase for output
    LS_amp = impro.crop_img(proper.prop_get_amplitude(wf), npupil, margin)\
            if get_amp is True else None
    LS_phase = impro.crop_img(proper.prop_get_phase(wf), npupil, margin)\
            if get_phase is True else None
    
    return wf, LS_amp, LS_phase
