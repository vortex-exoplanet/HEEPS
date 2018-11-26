from .apodization import apodization
from .vortex import vortex
from .lyotstop import lyotstop
from .lyot import lyot

def metis_modes(wfo, conf):
    """
    Various METIS modes Coronagraphic and non-coronagraphics:
        1. Vortex coronagraph (VC), to use put conf['MODE'] = 'VC'
        2. Ring apodized vortex coronagraph (RAVC), to use put conf['MODE'] = 'RAVC'
        3. Apodizing phase plate (APP), to use put conf['MODE'] = 'APP'
        4. No coronagraph, just Lyot-stop, to use put conf['MODE'] = 'OFFAXIS'    
        5. No coronagraph, but Ring apodizer and LS present, to use put conf['MODE'] = 'MASK'
        6. If conf['MODE'] = ELT, an ELT psf is generated    
        7. Classical Lyot coronagraph, to use put conf['MODE'] = 'CL' 
    """
    mode = conf['MODE']
    if mode == 'APP':
        lyotstop(wfo, conf, APP=True)
    elif mode == 'RAVC':
        RAVC = True
        apodization(wfo, conf, RAVC=RAVC)
        vortex(wfo, conf)
        lyotstop(wfo, conf, RAVC=RAVC)
    elif mode == 'CL': # classical lyot
        lyot(wfo, conf)
        lyotstop(wfo, conf)
    elif mode == 'VC':
        vortex(wfo, conf)
        lyotstop(wfo, conf)
    elif mode == 'OFFAXIS':
        print('No phase mask') # Lyot-stop (and ring apodizer) are present
        RAVC = False # False if VC, True if RAVC
        apodization(wfo, conf, RAVC=RAVC)
        lyotstop(wfo, conf, RAVC=RAVC)
    else:
        print('ELT PSF')
    return wfo
