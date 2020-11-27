from .circular_apodization import circular_apodization
import heeps.util.img_processing as impro
import proper
import numpy as np
import os.path
from astropy.io import fits 

def apodizer(wf, mode='RAVC', ravc_t=0.8, ravc_r=0.6, ravc_misalign=None, 
        ngrid=1024, npupil=285, file_ravc_amp='', file_ravc_phase='', 
        verbose=False, **conf):
    
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
    file_ravc_amp: str
    file_ravc_phase: str 
        ring apodizer files (optional)
    
    '''

    if mode in ['RAVC']:

        # load apodizer from files if provided
        if os.path.isfile(file_ravc_amp) and os.path.isfile(file_ravc_phase):
            if verbose is True:
                print('Load ring apodizer from files\n')
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
                print('Create ring apodizer')
                print('   ravc_t=%s, ravc_r=%s'\
                    %(round(ravc_t, 4), round(ravc_r, 4)))
                print('   ravc_misalign=%s'%np.around(ravc_misalign, 4))
                print('')

        # multiply the loaded apodizer
        proper.prop_multiply(wf, ring)
        
    return wf