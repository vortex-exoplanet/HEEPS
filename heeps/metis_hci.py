from .pupil import pupil
from .aberrations import wavefront_aberrations
from .coronagraphs import metis_modes
from .detector import detector

import numpy as np
from astropy.io import fits 
from multiprocessing import Pool
from copy import deepcopy
import time
import multiprocessing as mpro
from functools import partial
from sys import platform
import os

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
    wf = pupil(conf) 
    metis_modes(wf, conf)
    return 

def multi_cube(atm_screen,tip_tilt,conf,wf_start,iter):
    wf = deepcopy(wf_start)
    print(iter)
    if ((isinstance(atm_screen, (list, tuple, np.ndarray)) == True)):
        if (atm_screen.ndim == 3):
            atm_screen_iter = atm_screen[iter,:,:]
        else:
            atm_screen_iter = atm_screen
    if (tip_tilt.ndim == 2):
        tip_tilt_iter = tip_tilt[iter,:]
    else:
        tip_tilt_iter = tip_tilt
    conf['tip_tilt'] = tip_tilt_iter
    wavefront_aberrations(wf, AO_residuals=atm_screen_iter, **conf)
    metis_modes(wf, conf)
    psf = detector(wf, conf)
    return psf

def metis_hci(atm_screen, **conf):
    atm_screen = np.array(np.float32(atm_screen), ndmin=3)
    tip_tilt = np.array(conf['tip_tilt'])
    if (propagation_test(atm_screen,tip_tilt) == 'single'):
        wf = pupil(conf) 
        conf['tip_tilt'] = tip_tilt
        wavefront_aberrations(wf, AO_residuals=atm_screen, **conf)
        metis_modes(wf, conf)
        psf = detector(wf, conf)
        psf_cube = np.array(psf, ndmin=3)
    else:
        temp_prop(conf) # takes care of calibration files, if not present
        if (atm_screen.ndim == 3):
            length_cube = atm_screen.shape[0]
        if (tip_tilt.ndim == 2):
            length_cube = tip_tilt.shape[0]
        wf_start = pupil(conf)
        
        if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
            if conf['cpucount'] == None:
                conf['cpucount'] = mpro.cpu_count() - 1
            print('%s: simulation starts using %s cores.'\
                    %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), conf['cpucount']))
            p = Pool(conf['cpucount'])
            func = partial(multi_cube, atm_screen, tip_tilt, conf, wf_start)
            psf_cube = np.array(p.map(func, range(length_cube)))
        else:
            psf_cube = np.zeros((length_cube, conf['ndet'], conf['ndet']))
            for i in range(length_cube):
                psf_cube[i,:,:] = multi_cube(atm_screen,tip_tilt,conf,wf_start,i)
    
    return psf_cube
