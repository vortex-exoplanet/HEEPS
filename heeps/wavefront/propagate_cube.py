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
        vc_chrom_leak=2e-3, add_cl_det=False, add_cl_vort=False, tag=None, 
        onaxis=True, send_to=None, savefits=False, verbose=False, **conf):

    # update conf
    conf.update(mode=mode, ngrid=ngrid, cpu_count=cpu_count, 
        vc_chrom_leak=vc_chrom_leak, add_cl_det=add_cl_det, 
        add_cl_vort=add_cl_vort, tag=tag, onaxis=onaxis)

    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    if verbose == True:
        print('Create %s-axis PSF cube'%{True:'on',False:'off'}[onaxis])
        if add_cl_det is True:
            print('   adding chromatic leakage at detector plane: %s'%vc_chrom_leak)
        if add_cl_vort is True:
            print('   adding chromatic leakage at vortex plane: %s'%vc_chrom_leak)

    # preload amp screen if only one frame
    if len(amp_screens) == 1 and np.any(amp_screens) != None:
        proper.prop_multiply(wf1, pad_img(amp_screens, ngrid))
        # then create a cube of None values
        amp_screens = [None]*int(nframes/nstep + 0.5)
        if verbose == True:
            print('   preloading amplitude screen')

    # preload apodizer when no drift
    if ('APP' in mode) or ('RAVC' in mode and np.all(apo_misaligns) == None):
        wf1 = apodizer(wf1, verbose=False, **conf)
        conf['apo_loaded'] = True
        if verbose == True:
            print('   preloading %s apodizer'%mode)
    else:
        conf['apo_loaded'] = False

    # preload Lyot stop when no drift
    if ('VC' in mode or 'LC' in mode) and np.all(ls_misaligns) == None:
        conf['ls_mask'] = lyot_stop(wf1, apply_ls=False, verbose=False, **conf)
        if verbose == True:
            print('   preloading Lyot stop')

    # run simulation
    posvars = [phase_screens, amp_screens, tiptilts, apo_misaligns, ls_misaligns]
    kwargs = dict(verbose=False, **conf)
    psfs = multiCPU(propagate_one, posargs=[wf1], posvars=posvars, kwargs=kwargs,
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