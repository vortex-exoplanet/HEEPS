from .propagate_one import propagate_one
import multiprocessing as mpro
from functools import partial
from sys import platform
import numpy as np
import time
from astropy.io import fits 
import os.path

def propagate_cube(wf, phase_screens=None, tiptilts=None, misaligns=None, \
        case='', nframes=20, nstep=1, ndet=365, cpu_count=1, savefits=False, verbose=False, onaxis=True, **conf):

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
            print('   %s: e2e simulation starts, using %s cores.'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cpu_count))
        p = mpro.Pool(cpu_count)
        func = partial(propagate_one, wf, **conf)
        psfs = np.array(p.starmap(func, zip(phase_screens, tiptilts, misaligns)))
        p.close()
        p.join()
    else:
        print('      %s: e2e simulation starts, using 1 core.'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        psfs = np.zeros((nframes, ndet, ndet))
        for i, (phase_screen, tiptilt, misalign) \
                in enumerate(zip(phase_screens, tiptilts, misaligns)):
            psf = propagate_one(wf, conf, phase_screen, tiptilt, misalign)
            psfs[i,:,:] = psf                
    # if only one frame, make dim = 2
    if nframes == 1:
        psfs = psfs[0]
    # print elapsed time
    print('      elapsed %.3f seconds.'%(time.time() - t0))   
    print('')
    
    if savefits is True:
        prefix = '' if case == '' else '%s_'%case
        on_off = {True:'onaxis', False:'offaxis'}[onaxis]
        fits.writeto(os.path.join(conf['dir_output'], '%s%s_PSF_%s_%s.fits'\
            %(prefix, on_off, conf['band'], conf['mode'])), np.float32(psfs), overwrite=True)
    
    return psfs