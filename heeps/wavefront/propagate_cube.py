from .propagate_one import propagate_one
from heeps.optics import apodizer
from heeps.util.save2fits import save2fits
from heeps.util.notify import notify
from heeps.util.multiCPU import multiCPU
import numpy as np

def propagate_cube(wf, phase_screens, amp_screens, tiptilts, misaligns, cpu_count=1,
        vc_chrom_leak=2e-3, add_det_chrom_leak=False, add_vort_chrom_leak=False, 
        send_to=None, tag=None, onaxis=True, savefits=False, verbose=False, **conf):

    # preload amp screen if only one frame
    if len(amp_screens) == 1 and np.any(amp_screens) != None:
        import proper
        from heeps.util.img_processing import pad_img, resize_img
        amp_screens = np.nan_to_num(amp_screens[0])
        amp_screens = pad_img(resize_img(amp_screens, conf['npupil']), conf['ngrid'])
        proper.prop_multiply(wf, amp_screens)
        # then create a cube of None values
        amp_screens = [None]*int((conf['nframes']/conf['nstep']) + 0.5)

    # preload apodizer when no drift
    if np.all(misaligns) == None or 'APP' in conf['mode']:
        wf = apodizer(wf, verbose=False, **conf)
    
    if verbose == True:
        print('Create %s-axis PSF cube'%{True:'on',False:'off'}[onaxis])
        if add_det_chrom_leak is True:
            print('   adding chromatic leakage at detector plane: %s'%vc_chrom_leak)
        if add_vort_chrom_leak is True:
            print('   adding chromatic leakage at vortex plane: %s'%vc_chrom_leak)

    # run simulation
    posvars = [phase_screens, amp_screens, tiptilts, misaligns]
    kwargs = dict(onaxis=onaxis, verbose=False, **conf)
    psfs = multiCPU(propagate_one, posargs=[wf], posvars=posvars, kwargs=kwargs, \
        case='e2e simulation', cpu_count=cpu_count)

    # if only one wavefront, make dim = 2
    if len(psfs) == 1:
        psfs = psfs[0]

    # save cube of PSFs to fits file, and notify by email
    if savefits == True:
        tag = '' if tag is None else '%s_'%tag
        name = '%s%s_PSF'%(tag, {True: 'onaxis', False: 'offaxis'}[onaxis])
        filename = save2fits(psfs, name, **conf)
        notify('saved to %s'%filename, send_to)

    return psfs
