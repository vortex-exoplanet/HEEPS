from .propagate_one import propagate_one
from heeps.optics import apodizer, lyot_stop
from heeps.util.save2fits import save2fits
from heeps.util.notify import notify
from heeps.util.multiCPU import multiCPU
from heeps.util.img_processing import pad_img
from copy import deepcopy
import proper
import numpy as np

def propagate_cube(wf, phase_screens, amp_screens, tiptilts, apo_misaligns,
        ls_misaligns, nframes=10, nstep=1, mode='RAVC', ngrid=1024, cpu_count=1,
        vc_chrom_leak=2e-3, add_cl_det=False, add_cl_vort=False, avg=False, tag=None, 
        onaxis=True, send_to=None, savefits=False, verbose=False, **conf):

    # update conf
    conf.update(mode=mode, ngrid=ngrid, cpu_count=cpu_count, 
        vc_chrom_leak=vc_chrom_leak, add_cl_det=add_cl_det, 
        add_cl_vort=add_cl_vort, tag=tag, onaxis=onaxis)

    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    if verbose is True:
        print('Create %s-axis PSF cube'%{True:'on',False:'off'}[onaxis])
        if 'VC' in mode and onaxis is True:
            if add_cl_det is True:
                print('   adding chromatic leakage at detector plane: %s'%vc_chrom_leak)
            if add_cl_vort is True:
                print('   adding chromatic leakage at vortex plane: %s'%vc_chrom_leak)

    # preload amp screen if only one frame
    if np.any(amp_screens) != None:
        if len(amp_screens) == 1:
            proper.prop_multiply(wf1, pad_img(amp_screens[0], ngrid))
            # then create a cube of None values
            amp_screens = [None]*int(nframes/nstep + 0.5)
            if verbose == True:
                print('   preloading amplitude screen')
        elif verbose is True:
            print('   amplitude screen not preloaded: len(amp_screens)=%s'%len(amp_screens))

    # preload apodizer when no drift
    conf['apo_loaded'] = False
    if ('APP' in mode) or ('RAVC' in mode):
        if np.all(apo_misaligns[:-1] == apo_misaligns[1:]):
            wf1 = apodizer(wf1, verbose=verbose, **conf)
            conf['apo_loaded'] = True
        elif verbose is True:
            print('   %s apodizer not preloaded: len(apo_misalign)=%s'%(mode, len(apo_misaligns)))

    # preload Lyot stop when no drift
    if mode in ['CVC', 'RAVC', 'CLC', 'IMG', 'LMS']:
        if np.all(ls_misaligns[:-1] == ls_misaligns[1:]):
            conf['ls_mask'] = lyot_stop(wf1, apply_ls=False, verbose=verbose, **conf)
        elif verbose is True:
            print('   Lyot stop not preloaded: len(ls_misalign)=%s'%len(ls_misaligns))

    # run simulation
    del conf['apo_misalign'], conf['ls_misalign']
    posvars = [phase_screens, amp_screens, tiptilts, apo_misaligns, ls_misaligns]
    kwargs = dict(verbose=False, **conf)
    psfs = multiCPU(propagate_one, posargs=[wf1], posvars=posvars, kwargs=kwargs,
        case='e2e simulation', cpu_count=cpu_count)

    # optional: average cube
    if avg is True:
        if verbose is True:
            print('Average PSF cube\n')
        psfs = np.mean(psfs, axis=0)
    elif verbose is True:
        print('')

    # save cube of PSFs to fits file, and notify by email
    if savefits == True:
        tag = '' if tag is None else '%s_'%tag
        name = '%s%s_PSF'%(tag, {True: 'onaxis', False: 'offaxis'}[onaxis])
        filename = save2fits(psfs, name, **conf)
        notify('saved to %s'%filename, send_to)

    return psfs