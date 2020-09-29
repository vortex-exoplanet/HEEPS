from .propagate_one import propagate_one
from heeps.util.save2fits import save2fits
from heeps.util.notify import notify
import multiprocessing as mpro
from functools import partial
from sys import platform
import numpy as np
import time

def propagate_cube(wf, phase_screens=None, tiptilts=None, misaligns=None,
        dir_output='output_files', case='', band='L', mode='RAVC', nframes=20, 
        nstep=1, ndet=365, cpu_count=1, savefits=False, send_to=None, 
        verbose=False, onaxis=True, **conf):

    # update conf
    conf.update(ndet=ndet)
    
    # calculate number of wavefronts
    nwf = int((nframes/nstep) + 0.5)
    phase_screens = phase_screens[:nframes][::nstep] if np.any(phase_screens) else [None]*nwf
    tiptilts = tiptilts[:nframes][::nstep] if np.any(tiptilts) else [None]*nwf
    misaligns = misaligns[:nframes][::nstep] if np.any(misaligns) else [None]*nwf
    
    # run simulation
    t0 = time.time()
    if cpu_count != 1 and platform in ['linux', 'linux2', 'darwin']:
        if cpu_count == None:
            cpu_count = mpro.cpu_count() - 1
        if verbose is True:
            print('%s: e2e simulation starts, using %s cores'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cpu_count))
        p = mpro.Pool(cpu_count)
        func = partial(propagate_one, wf, **conf)
        psfs = np.array(p.starmap(func, zip(phase_screens, tiptilts, misaligns)))
        p.close()
        p.join()
    else:
        if verbose is True:
            print('%s: e2e simulation starts, using 1 core'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        psfs = np.zeros((nwf, ndet, ndet))
        for i, (phase_screen, tiptilt, misalign) \
                in enumerate(zip(phase_screens, tiptilts, misaligns)):
            psf = propagate_one(wf, phase_screen, tiptilt, misalign, **conf)
            psfs[i,:,:] = psf                
    if verbose is True:
        print('%s: finished, elapsed %.3f seconds\n'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), time.time() - t0))

    # if only one wavefront, make dim = 2
    if nwf == 1:
        psfs = psfs[0]

    # save cube of PSFs to fits file, and notify by email
    if savefits == True:
        prefix = '' if case == '' else '%s_'%case
        name = '%s%s_PSF'%(prefix, {True: 'onaxis', False: 'offaxis'}[onaxis])
        filename = save2fits(psfs, name, **conf)
        notify('%s  created.'%filename, send_to)

    return psfs
