from heeps.optics import apodizer, fp_mask, lyot_stop, detector
from heeps.util.img_processing import resize_img, pad_img
from copy import deepcopy
import numpy as np
import proper
from astropy.io import fits
import os.path

def propagate_one(wf, phase_screen=None, tiptilt=None, misalign=None, 
        npupil=285, ngrid=1024, savefits=False, onaxis=True, verbose=False, **conf):
            
    """ 
    Propagate one single wavefront.
    An off-axis PSF can be obtained by switching onaxis to False,
    thereby decentering the focal plane mask (if any).
    """
    # update conf
    conf.update(
        ravc_misalign=misalign,
        ngrid=ngrid)

    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    # apply phase screen (scao residuals, ncpa, petal piston)
    if phase_screen is not None:
        phase_screen = np.nan_to_num(phase_screen)# in microns
        phase_screen = resize_img(phase_screen, npupil)
        phase_screen = pad_img(phase_screen, ngrid)        
        proper.prop_add_phase(wf1, phase_screen)

    # apply tip-tilt (Zernike 2,3)
    if tiptilt is not None:
        proper.prop_zernikes(wf1, [2,3], tiptilt)

    # pupil-plane apodization
    wf1, apo_amp, apo_phase = apodizer(wf1, get_amp=True, verbose=verbose, **conf)
    # focal-plane mask, only in 'on-axis' configuration
    if onaxis == True:
        wf1 = fp_mask(wf1, verbose=verbose, **conf)
    # Lyot-stop or APP
    wf1, ls_amp, ls_phase = lyot_stop(wf1, get_amp=True, verbose=verbose, **conf)
    # detector
    psf = detector(wf1, verbose=verbose, **conf)

    # save psf as fits file
    if savefits == True:
        on_off = {True: 'onaxis', False: 'offaxis'}
        fits.writeto(os.path.join(conf['dir_output'], '%s_PSF_%s_%s.fits'\
            %(on_off[onaxis], conf['band'], conf['mode'])), np.float32(psf), overwrite=True)

    return psf
