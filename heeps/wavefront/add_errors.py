from heeps.util.img_processing import pad_img
import proper
import numpy as np

def add_errors(wf, phase_screen=None, amp_screen=None, 
        tiptilt=None, ngrid=1024, zern=None, **conf):

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

    # add constant zernike
    if zern is not None:
        proper.prop_zernikes(wf, [zern[0]], [zern[1]])

    return wf