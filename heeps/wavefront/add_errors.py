from heeps.optics import apodizer
from heeps.util.img_processing import pad_img
import proper
import numpy as np

def add_errors(wf, onaxis=True, phase_screen=None, amp_screen=None, 
        tiptilt=None, apo_misalign=None, ngrid=1024,
        apo_loaded=False, verbose=False, **conf):

    # apply phase screen (scao residuals, ncpa, petal piston)
    if phase_screen is not None:
        assert phase_screen.ndim == 2, "phase_screen dim must be 2."
        proper.prop_add_phase(wf, pad_img(phase_screen, ngrid))

    # apply amplitude screen (Talbot effect)
    if amp_screen is not None:
        assert amp_screen.ndim == 2, "amp_screen dim must be 2."
        proper.prop_multiply(wf, pad_img(amp_screen, ngrid))

    # apply tip-tilt (Zernike 2,3)
    if tiptilt is not None:
        proper.prop_zernikes(wf, [2,3], np.array(tiptilt, ndmin=1))

    # pupil-plane apodization (already preloaded if no RA misalign)
    if apo_loaded == False:
        wf = apodizer(wf, ngrid=ngrid, apo_misalign=apo_misalign,
            onaxis=onaxis, verbose=verbose, **conf)

    return wf