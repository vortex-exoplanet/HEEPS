from .propagate_one import propagate_one
from heeps.optics import apodizer
from heeps.util.save2fits import save2fits
from heeps.util.notify import notify
import multiprocessing as mpro
from functools import partial
from sys import platform
import numpy as np
import time

def propagate_cube(wf, phase_screens, amp_screens, tiptilts, misaligns,
        cpu_count=1, send_to=None, tag=None, onaxis=True, savefits=False, 
        verbose=False, **conf):

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
        print('Create cube of %s-axis PSFs'%{True:'on',False:'off'}[onaxis])

    # run simulation
    t0 = time.time()
    if cpu_count != 1 and platform in ['linux', 'linux2', 'darwin']:
        if cpu_count == None:
            cpu_count = mpro.cpu_count() - 1
        if verbose is True:
            print('   %s: e2e simulation starts, using %s cores'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cpu_count))
        p = mpro.Pool(cpu_count)
        func = partial(propagate_one, wf, onaxis=onaxis, verbose=False, **conf)
        psfs = np.array(p.starmap(func, zip(phase_screens, amp_screens, tiptilts, misaligns)))
        p.close()
        p.join()
    else:
        if verbose is True:
            print('   %s: e2e simulation starts, using 1 core'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        for i, (phase_screen, amp_screen, tiptilt, misalign) \
                in enumerate(zip(phase_screens, amp_screens, tiptilts, misaligns)):
            psf = propagate_one(wf, phase_screen, amp_screen, tiptilt, misalign, \
                onaxis=onaxis, verbose=False, **conf)
            psfs = psf if i == 0 else np.dstack((psfs.T, psf.T)).T
    if verbose is True:
        print('   %s: finished, elapsed %.3f seconds'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), time.time() - t0))

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
