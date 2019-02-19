from .apodization import apodization
from .vortex import vortex
from .lyotstop import lyotstop
from .lyot import lyot

def metis_modes(wf, conf):
    """ Various METIS HCI modes:
    0. conf['mode'] = 'ELT'  for no coronagraph (only telescope)
    1. conf['mode'] = 'RAVC' for Ring Apodized Vortex Coronagraph
    2. conf['mode'] = 'CVC'  for Classical Vortex Coronagraph
    3. conf['mode'] = 'APP'  for Apodizing Phase Plate
    4. conf['mode'] = 'CLC'  for Classical Lyot Coronagraph
    
    An off-axis PSF can be obtained by switching conf['onaxis'] to False,
    thereby decentering the focal plane mask (if any).
    """
    
    if conf['mode'] == 'ELT': 
        conf['LS_params'] = [1., -0.3, 0.]  # means no Lyot stop
        lyotstop(wf, conf)
    elif conf['mode'] == 'RAVC':
        RAVC = True                         # ring apodizer
        apodization(wf, conf, RAVC=RAVC)
        if conf['onaxis'] == True:
            vortex(wf, conf)
        lyotstop(wf, conf, RAVC=RAVC)
    elif conf['mode'] == 'CVC':
        if conf['onaxis'] == True:
            vortex(wf, conf)
        lyotstop(wf, conf)
    elif conf['mode'] == 'APP':
        APP = True                          # apodizing phase plate
        if conf['onaxis'] == True:
            lyotstop(wf, conf, APP=APP)
    elif conf['mode'] == 'CLC':
        conf['LS_params'] = [0.8, 0.1, 1.1] # CLC lyot-stop parameters
        if conf['onaxis'] == True:
            lyot(wf, conf)
        lyotstop(wf, conf)
    
    return wf
