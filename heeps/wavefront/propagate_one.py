from heeps.optics import apodizer, fp_mask, lyot_stop, detector
from heeps.util.img_processing import resize_img, pad_img
from heeps.util.save2fits import save2fits
from copy import deepcopy
import proper
import numpy as np

def propagate_one(wf, phase_screen=None, amp_screen=None, tiptilt=None, misalign=[0,0,0,0,0,0], 
        ngrid=1024, npupil=285, tag=None, onaxis=True, savefits=False, verbose=False, 
        **conf):
            
    """ 
    Propagate one single wavefront.
    An off-axis PSF can be obtained by switching onaxis to False,
    thereby decentering the focal plane mask (if any).
    """

    # update conf
    conf.update(ngrid=ngrid, npupil=npupil)
    
    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    # apply phase screen (scao residuals, ncpa, petal piston)
    if phase_screen is not None:
        assert phase_screen.ndim == 2, "phase_screen dim must be 2."
        phase_screen = np.nan_to_num(phase_screen)
        phase_screen = pad_img(resize_img(phase_screen, npupil), ngrid)
        proper.prop_add_phase(wf1, phase_screen)

    # apply amplitude screen (Talbot effect)
    if amp_screen is not None:
        assert amp_screen.ndim == 2, "amp_screen dim must be 2."
        amp_screen = np.nan_to_num(amp_screen)
        amp_screen = pad_img(resize_img(amp_screen, npupil), ngrid)
        proper.prop_multiply(wf1, amp_screen)

    # apply tip-tilt (Zernike 2,3)
    if tiptilt is not None:
        # # translate the tip/tilt from lambda/D into RMS phase errors
        # # RMS = x/2 = ((lam/D)*(D/2))/2 = lam/4
        # tiptilt = np.array(tiptilt, ndmin=1)*conf['lam']/4
        # translate the tip/tilt from mas into RMS phase errors
        import astropy.units as u
        tiptilt = np.array(tiptilt, ndmin=1)*u.mas.to('rad')*conf['diam_ext']/4
        proper.prop_zernikes(wf1, [2,3], tiptilt)
    
    # pupil-plane apodization: if misalign set to None, apodizer already preloaded
    if misalign is not None:
        conf.update(ravc_misalign=misalign)
        wf1 = apodizer(wf1, verbose=verbose, **conf)
    
    # focal-plane mask, only in 'on-axis' configuration
    if onaxis == True:
        wf1 = fp_mask(wf1, verbose=verbose, **conf)
    
    # Lyot-stop
    wf1 = lyot_stop(wf1, verbose=verbose, **conf)
    
    # detector
    psf = detector(wf1, verbose=verbose, **conf)

    # save psf as fits file
    if savefits == True:
        tag = '' if tag is None else '%s_'%tag
        name = '%s%s_PSF'%(tag, {True: 'onaxis', False: 'offaxis'}[onaxis])
        save2fits(psf, name, **conf)

    return psf
