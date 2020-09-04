from .propagate_one import propagate_one
import multiprocessing as mpro
from functools import partial
from sys import platform
import numpy as np
import time
from astropy.io import fits 
import os.path

def propagate_cube(wf, conf, phase_screens=None, misaligns=None, zernikes=None, \
        case='', savefits=False, verboses=False):

    nframes = conf['nframes']
    nstep = conf['nstep']
    nwf = int((nframes/nstep) + 0.5)
    phase_screens = phase_screens[:nframes][::nstep] if np.any(phase_screens) else [None]*nwf
    misaligns = misaligns[:nframes][::nstep] if np.any(misaligns) else [None]*nwf
    zernikes = zernikes[:nframes][::nstep] if np.any(zernikes) else [None]*nwf

    t0 = time.time()
    if conf['cpu_count'] != 1 and platform in ['linux', 'linux2', 'darwin']:
        if conf['cpu_count'] == None:
            conf['cpu_count'] = mpro.cpu_count() - 1
        print('      %s: e2e simulation starts, using %s cores.'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), conf['cpu_count']))
        p = mpro.Pool(conf['cpu_count'])
        func = partial(propagate_one, wf, conf)
        psfs = np.array(p.starmap(func, zip(phase_screens, misaligns, zernikes)))
        p.close()
        p.join()
    else:
        print('      %s: e2e simulation starts, using 1 core.'\
                %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        psfs = np.zeros((conf['nframes'], conf['ndet'], conf['ndet']))
        for i, (phase_screen, misalign, \
                zernike) in enumerate(zip(phase_screens, misaligns, zernikes)):
            psf = propagate_one(wf, conf, phase_screen, misalign, zernike)
            psfs[i,:,:] = psf                
    # if only one frame, make dim = 2
    if conf['nframes'] == 1:
        psfs = psfs[0]
    # print elapsed time
    print('      elapsed %.3f seconds.'%(time.time() - t0))   
    print('')
    
    if savefits is True:
        conf['prefix'] = '%s_'%case
        on_off = {True:'onaxis', False:'offaxis'}[conf['onaxis']]
        filename = '%s%s'%(conf['prefix'], on_off)+'_%s_'+'%s_%s'%(conf['band'], conf['mode'])
        fits.writeto(os.path.join(conf['dir_output'], filename%'PSF') \
                + '.fits', np.float32(psfs), overwrite=True)
    
    return psfs