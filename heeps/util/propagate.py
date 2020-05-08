import numpy as np
from copy import deepcopy
from heeps.aberrations import wavefront_aberrations
from heeps.coronagraphs import apodization, vortex, lyotstop, lyot
from heeps.detector import detector

def prop_cube(wf_start, conf, atm_screens, ncpa_screens, 
            piston_screens, misaligns, zernikes):
    pass

def prop_one(wf_start, conf, atm_screen, ncpa_screen, 
            piston_screen, misalign, zernike):
    """ Create a function to propagate one single wavefront """

    # keep a copy of the input pupil wavefront
    wf = deepcopy(wf_start)
    # apply wavefront aberrations
    wavefront_aberrations(wf, zernike=zernike, piston_screen=piston_screen, \
            atm_screen=atm_screen, ncpa_screen=ncpa_screen, **conf) #resized=True, nans=False, **conf)
    conf['RAVC_misalign'] = misalign
    # define flags for special cases
    RAVC = True if conf['mode'] in ['RAVC'] else False
    APP = True if conf['mode'] in ['APP'] and conf['onaxis'] is True else False 
    # pupil-plane apodization
    wf, apo_amp, apo_phase = apodization(wf, conf, RAVC=RAVC)
    # select focal-plane mask, and Lyot-stop
    if conf['mode'] in ['ELT', 'APP']:
        conf['LS_params'] = [1., -0.313, 0.] # no Lyot stop
    elif conf['mode'] in ['CLC']:
        conf['LS_params'] = [0.8, 0.1, 1.1]
        if conf['onaxis'] == True:
            lyot(wf, conf)
    elif conf['mode'] in ['CVC', 'RAVC']:
        conf['LS_params'] = [0.98, 0.03, 1.1]
        if conf['onaxis'] == True:
            vortex(wf, conf)
    # lyot-stop
    wf, LS_amp, LS_phase = lyotstop(wf, conf, RAVC=RAVC, APP=APP)
    # get science image
    conf['ndet'] = int(np.ceil(2*conf['hfov']*1000/conf['pscale']))    
    if conf['ndet'] % 2:
        conf['ndet'] += 1
    psf = detector(wf, conf)
    # release memory: global variable set to None
    wf = None
    
    if conf['full_output'] is True:
        return psf, apo_amp, LS_amp
    else:
        return psf