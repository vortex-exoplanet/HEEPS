from heeps.optics import apodizer, fp_mask, lyot_stop, detector
from heeps.wavefront.add_errors import add_errors
from heeps.util.save2fits import save2fits
from copy import deepcopy
import proper
import numpy as np

def propagate_one(wf, phase_screen=None, amp_screen=None, tiptilt=None,
        apo_misalign=None, ls_misalign=None, mode='RAVC', apo_loaded=False,
        ngrid=1024, npupil=285, vc_chrom_leak=2e-3, add_cl_det=False, tag=None,
        fp_offsets=None, onaxis=True, savefits=False, verbose=False, **conf):

    """
    Propagate one single wavefront.
    An off-axis PSF can be obtained by switching onaxis to False,
    thereby decentering the focal plane mask (if any).
    """

    if verbose == True:
        print('Create single %s-axis PSF'%{True:'on',False:'off'}[onaxis])

    # update conf
    conf.update(mode=mode, ngrid=ngrid, npupil=npupil, tag=tag, onaxis=onaxis,
            vc_chrom_leak=vc_chrom_leak, add_cl_det=add_cl_det)
    # no misalignment for off-axis PSF
    if onaxis == True:
        conf.update(apo_misalign=apo_misalign, ls_misalign=ls_misalign)

    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    # pupil-plane apodization (RAP, APP) already preloaded if no RA drift
    if apo_loaded == False:
        wf = apodizer(wf, verbose=verbose, **conf)

    # apply wavefront errors (SCAO residuals, NCPA, Talbot effect, ...)
    wf1 = add_errors(wf1, phase_screen=phase_screen, amp_screen=amp_screen,
        tiptilt=tiptilt, verbose=verbose, **conf)

    # imaging a point source
    def point_source(wf1, verbose, conf):
        if onaxis == True:
            if add_cl_det is True and 'VC' in mode: # add chromatic leakage (vortex only)
                cl = deepcopy(wf1)
                cl._wfarr = np.flip(cl._wfarr) # 2 FFTs
                cl = lyot_stop(cl, verbose=verbose, **conf)
                chrom_leak = cl._wfarr*np.sqrt(vc_chrom_leak)
            else:
                chrom_leak = 0
            wf1 = fp_mask(wf1, verbose=verbose, **conf) # focal-plane mask (onaxis only)
            wf1 = lyot_stop(wf1, verbose=verbose, **conf)
            wf1._wfarr += chrom_leak
        else:
            wf1._wfarr = np.flip(wf1._wfarr) # 2 FFTs
            wf1 = lyot_stop(wf1, verbose=verbose, **conf)
        return detector(wf1, verbose=verbose, **conf)

    # imaging a point source
    if fp_offsets is None:
        psf = point_source(wf1, verbose, conf)
    # imaging a finite size star
    else:
        psf = 0
        for offset in fp_offsets:
            point = deepcopy(wf1)
            proper.prop_zernikes(point, [2,3], np.array(offset, ndmin=1))
            psf += point_source(point, False, conf)
        psf /= len(fp_offsets)

    # save psf as fits file
    if savefits == True:
        tag = '' if tag is None else '%s_'%tag
        name = '%s%s_PSF'%(tag, {True: 'onaxis', False: 'offaxis'}[onaxis])
        save2fits(psf, name, **conf)

    return psf