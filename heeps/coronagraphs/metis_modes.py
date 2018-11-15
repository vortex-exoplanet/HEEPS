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
    conf['PHASE_APODIZER_FILE'] = 0

    if mode == 'APP': 
        conf['PHASE_APODIZER_FILE'] = 'app_phase_cut.fits' 
        lyotstop(wfo, conf, RAVC=False)   
    if mode == 'RAVC':    
        apodization(wfo, conf, RAVC=True)
        vortex(wfo, conf)
        lyotstop(wfo, conf, RAVC=True)
    elif mode == 'CL': # classical lyot
        lyot(wfo, conf)
        lyotstop(wfo, conf, RAVC=False)
    elif mode == 'VC':
        vortex(wfo, conf)
        lyotstop(wfo, conf, RAVC=False)
    elif mode == 'OFFAXIS':
        print('No Coronagraph')    
        lyotstop(wfo, conf, RAVC=False)
    elif mode == 'MASK':
        print('Ring apodizer and LS present')    
        apodization(wfo, conf, RAVC=True)
        lyotstop(wfo, conf, RAVC=True)
    else: 
        print('ELT PSF')
    return wfo
