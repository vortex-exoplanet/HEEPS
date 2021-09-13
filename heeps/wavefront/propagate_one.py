from heeps.optics import fp_mask, lyot_stop, detector
from heeps.wavefront.add_errors import add_errors
from heeps.util.save2fits import save2fits
from copy import deepcopy
import proper
import numpy as np

def propagate_one(wf, phase_screen=None, amp_screen=None, tiptilt=None, misalign=[0,0,0,0,0,0], 
        ngrid=1024, npupil=285, vc_chrom_leak=2e-3, add_cl_det=False, fp_offsets=None, 
        tag=None, onaxis=True, savefits=False, verbose=False, **conf):
            
    """ 
    Propagate one single wavefront.
    An off-axis PSF can be obtained by switching onaxis to False,
    thereby decentering the focal plane mask (if any).
    """

    if verbose == True:
        print('Create single %s-axis PSF'%{True:'on',False:'off'}[onaxis])

    # update conf
    conf.update(ngrid=ngrid, npupil=npupil, vc_chrom_leak=vc_chrom_leak, \
            add_cl_det=add_cl_det, tag=tag, onaxis=onaxis)
    
    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    # apply phase screen (scao residuals, ncpa, petal piston)
    wf1 = add_errors(wf1, phase_screen=phase_screen, amp_screen=amp_screen, \
        tiptilt=tiptilt, misalign=misalign, verbose=verbose, **conf)
    
    # imaging a point source
    def point_source(wf1, verbose, conf):
        if onaxis == True: # focal-plane mask, only in 'on-axis' configuration
            if add_cl_det is True:
                cl = deepcopy(wf1)
                cl._wfarr = np.flip(cl._wfarr) # 2 FFTs
                cl = lyot_stop(cl, verbose=verbose, **conf)
            wf1 = fp_mask(wf1, verbose=verbose, **conf)
            wf1 = lyot_stop(wf1, verbose=verbose, **conf)
            if add_cl_det is True:
                wf1._wfarr += cl._wfarr*np.sqrt(vc_chrom_leak)
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
