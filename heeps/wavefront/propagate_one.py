from heeps.optics import apodizer, fp_mask, lyot_stop, detector
from heeps.util.img_processing import pad_img
from heeps.util.save2fits import save2fits
from copy import deepcopy
import proper
import numpy as np

def propagate_one(wf, phase_screen=None, amp_screen=None, tiptilt=None, misalign=[0,0,0,0,0,0], 
        ngrid=1024, npupil=285, fp_offsets=None, tag=None, onaxis=True, savefits=False, 
        verbose=False, **conf):
            
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
        proper.prop_add_phase(wf1, pad_img(phase_screen, ngrid))

    # apply amplitude screen (Talbot effect)
    if amp_screen is not None:
        assert amp_screen.ndim == 2, "amp_screen dim must be 2."
        proper.prop_multiply(wf1, pad_img(amp_screen, ngrid))

    # apply tip-tilt (Zernike 2,3)
    if tiptilt is not None:
        proper.prop_zernikes(wf1, [2,3], np.array(tiptilt, ndmin=1))
    
    if verbose == True:
        print('Create single %s-axis PSF'%{True:'on',False:'off'}[onaxis])

    # pupil-plane apodization: if RA misalign set to None (or APP cube calc),
    # then apodizer was already preloaded
    if misalign is not None or ('APP' in conf['mode'] and onaxis == False):
        conf.update(ravc_misalign=misalign)
        wf1 = apodizer(wf1, onaxis=onaxis, verbose=verbose, **conf)

    # imaging a point source
    def point_source(wfo, verbose, conf):
        if onaxis == True: # focal-plane mask, only in 'on-axis' configuration
            wfo = fp_mask(wfo, verbose=verbose, **conf)
        wfo = lyot_stop(wfo, verbose=verbose, **conf)
        return detector(wfo, verbose=verbose, **conf)
    if fp_offsets is None:
        psf = point_source(wf1, verbose, conf)
    
    # imaging a finite size star
    else:
        psf = 0
        for i,offset in enumerate(fp_offsets):
            wfo = deepcopy(wf1)
            proper.prop_zernikes(wfo, [2,3], np.array(offset, ndmin=1))
            psf += point_source(wfo, False, conf)
        psf /= len(fp_offsets)

    # save psf as fits file
    if savefits == True:
        tag = '' if tag is None else '%s_'%tag
        name = '%s%s_PSF'%(tag, {True: 'onaxis', False: 'offaxis'}[onaxis])
        save2fits(psf, name, **conf)

    return psf
