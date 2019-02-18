from .pupil import pupil
from .aberrations import wavefront_aberrations
from .coronagraphs import metis_modes
from .detector import detector


from copy import deepcopy
import numpy as np
from astropy.io import fits 
from multiprocessing import Pool
from functools import partial
from sys import platform

def check_dim(input):
    r = 0
    try:
        r = input.ndim
    except: AttributeError
    return r

def propagation_test(a,b):
    if check_dim(a)==3 or check_dim(b)==2:
        out = 'multi'
    else:
        out = 'single'
    return out

def temp_prop(conf):
    wfo = pupil(conf) 
    metis_modes(wfo, conf)
    return 

def multi_cube(atm_screen,TILT,conf,wfo1,iter):
    wfo = deepcopy(wfo1)  
    print(iter)
    if ((isinstance(atm_screen, (list, tuple, np.ndarray)) == True)):
        if (atm_screen.ndim == 3):
            atm_screen_iter = atm_screen[iter,:,:]
        else:
            atm_screen_iter = atm_screen
    if (TILT.ndim == 2):
        TILT_iter = TILT[iter,:]
    else:
        TILT_iter = TILT
    wavefront_aberrations(wfo, AO_residuals=atm_screen_iter, tip_tilt=TILT_iter, **conf)
    metis_modes(wfo, conf)
    psf = detector(wfo, conf)
    return psf

def metis_hci(conf, atm_screen, TILT):
    if (propagation_test(atm_screen,TILT) == 'single'):
        wfo = pupil(conf) 
        wavefront_aberrations(wfo,  AO_residuals=atm_screen, tip_tilt=TILT, **conf)
        metis_modes(wfo, conf)
        psf = detector(wfo, conf)
        fits.writeto(conf['OUT_DIR'] + conf['PREFIX'] + '_PSF_'+ conf['MODE'] +'.fits', psf, overwrite=True)
    else:
        temp_prop(conf) # takes care of calibration files, if not present
        if (atm_screen.ndim == 3):
            length_cube = atm_screen.shape[0]
        if (TILT.ndim == 2):
            length_cube = TILT.shape[0]
        wfo1 = pupil(conf)
        
        if conf['cpucount'] != 1 and platform in ('linux', 'linux2', 'darwin'):
            p = Pool(conf['cpucount'])
            func = partial(multi_cube, atm_screen, TILT, conf, wfo1)
            psf_cube = np.array(p.map(func, range(length_cube)))
        else:
            psf_cube = np.zeros((length_cube, conf['N_D'], conf['N_D']))
            for i in range(length_cube):
                psf_cube[i,:,:] = multi_cube(atm_screen,TILT,conf,wfo1,i)
        fits.writeto(conf['OUT_DIR'] + conf['PREFIX'] + '_PSF_cube_'+ conf['MODE'] +'.fits', psf_cube, overwrite=True)
        psf = psf_cube[0]
    
    return psf


