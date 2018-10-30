from .apodization import apodization
from .vortex import vortex
from .lyotstop import lyotstop
from .lyot import lyot

def metis_modes(wfo, conf):
    mode = conf['MODE']
    if mode == 'APP': 
        RAVC = False    
        lyotstop(wfo, conf, RAVC)
    
    conf['PHASE_APODIZER_FILE'] = 0
    if mode == 'RAVC':    
        RAVC = True
        apodization(wfo, conf, RAVC=True)
        vortex(wfo, conf)
        lyotstop(wfo, conf, RAVC)
    elif conf['MODE'] == 'CL': # classical lyot
        lyot(wfo, conf)
        _, LS_pupil = lyotstop(wfo, conf, RAVC=False)
    elif mode == 'VC':
        RAVC = False
        vortex(wfo, conf)
        lyotstop(wfo, conf, RAVC)
    elif mode == 'OFFAXIS':
        print('No Coronagraph')    
        RAVC = False
        lyotstop(wfo, conf, RAVC)
    elif mode == 'MASK':
        print('Ring apodizer and LS present')    
        RAVC = True
        apodization(wfo, conf, RAVC=True)
        lyotstop(wfo, conf, RAVC)
    else:
        print('ELT PSF')    
    return wfo






