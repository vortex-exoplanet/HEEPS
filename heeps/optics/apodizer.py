from .circular_apodization import circular_apodization
import heeps.util.img_processing as impro
import proper
import numpy as np
import os.path
from astropy.io import fits 

def apodizer(wf, mode='RAVC', ravc_t=0.8, ravc_r=0.6, ravc_misalign=None, 
        ngrid=1024, npupil=285, file_app_phase='', file_app_amp='', 
        file_ravc_amp='', file_ravc_phase='', onaxis=True, verbose=False, **conf):
    
    ''' Create a wavefront object at the entrance pupil plane. 
    The pupil is either loaded from a fits file, or created using 
    pupil parameters.
    Can also select only one petal and mask the others.

    wf: WaveFront
        PROPER wavefront object
    mode: str
        HCI mode
    ravc_t: float
        RA transmittance
    ravc_r: float
        RA radius
    ravc_misalign: list of float
        RA misalignment
    ngrid: int
        number of pixels of the wavefront array
    npupil: int
        number of pixels of the pupil
    file_app_amp: str
    file_app_phase: str 
        apodizing phase plate files
    file_ravc_amp: str
    file_ravc_phase: str 
        ring apodizer files (optional)

    '''

    # case 1: Ring Apodizer
    if 'RAVC' in mode and ravc_r > 0:

        # load apodizer from files if provided
        if os.path.isfile(file_ravc_amp) and os.path.isfile(file_ravc_phase):
            if verbose is True:
                print('   apply ring apodizer from files')
            # get amplitude and phase data
            RAVC_amp = fits.getdata(file_ravc_amp)
            RAVC_phase = fits.getdata(file_ravc_phase)
            # resize to npupil
            RAVC_amp = impro.resize_img(RAVC_amp, npupil)
            RAVC_phase = impro.resize_img(RAVC_phase, npupil)
            # pad with zeros to match PROPER gridsize
            RAVC_amp = impro.pad_img(RAVC_amp, ngrid)
            RAVC_phase = impro.pad_img(RAVC_phase, ngrid)
            # build complex apodizer
            apo = RAVC_amp*np.exp(1j*RAVC_phase)

        # or else, define the apodizer as a ring (with % misalignments)
        else:
            # RAVC misalignments
            ravc_misalign = [0,0,0,0,0,0] if ravc_misalign is None else list(ravc_misalign)
            dx_amp, dy_amp, dz_amp = ravc_misalign[0:3]
            dx_phase, dy_phase, dz_phase = ravc_misalign[3:6]
            # create apodizer
            ring = circular_apodization(wf, ravc_r, 1., ravc_t, xc=dx_amp, \
                yc=dy_amp, NORM=True)
            if verbose is True:
                print('   apply ring apodizer: ravc_t=%s, ravc_r=%s'\
                    %(round(ravc_t, 4), round(ravc_r, 4)))

        # multiply the loaded apodizer
        proper.prop_multiply(wf, ring)


    # case 2: Apodizing Phase Plate
    elif 'APP' in mode:
        # get amplitude and phase data
        if os.path.isfile(file_app_amp):
            if verbose is True:
                print('   apply APP stop (amplitude)')
            APP_amp = fits.getdata(file_app_amp)
        else:
            APP_amp = np.ones((npupil, npupil))
        if os.path.isfile(file_app_phase) and onaxis == True:
            if verbose is True:
                print('   apply APP phase')
            APP_phase = fits.getdata(file_app_phase)
        else:
            APP_phase = np.zeros((npupil, npupil))
        # resize to npupil
        APP_amp = impro.resize_img(APP_amp, npupil)
        APP_phase = impro.resize_img(APP_phase, npupil)
        # rotate for negative PSF
        if 'neg' in mode:
            APP_amp = np.rot90(APP_amp, 2)
            APP_phase = np.rot90(APP_phase, 2)
        # pad with zeros to match PROPER ngrid
        APP_amp = impro.pad_img(APP_amp, ngrid, 0)
        APP_phase = impro.pad_img(APP_phase, ngrid, 0)
        
        # multiply the loaded APP
        proper.prop_multiply(wf, APP_amp*np.exp(1j*APP_phase))
            
    return wf