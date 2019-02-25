import proper
import numpy as np
from astropy.io import fits
from .atmosphere import atmosphere
from .island_effect_piston import island_effect_piston
from .static_ncpa import static_ncpa

def wavefront_aberrations(wf, tip_tilt=None, petal_piston=None, 
            atm_screen=None, static_ncpa=False, polish_error=False, **conf):
    
    # add zernikes 2: X tilt, and 3: Y tilt
    if tip_tilt is not None:
        # translate the tip/tilt from lambda/D into RMS phase errors
        lam = conf['lam']
        tip_tilt = np.array(tip_tilt)*lam/4
        proper.prop_zernikes(wf, [2,3], tip_tilt)
    
    # add petal piston (island effect)
#    if petal_piston is not None:
#        wf = island_effect_piston(wf, petal_piston, **conf)
    
    # add atmosphere
    if atm_screen is not None:
        wf = atmosphere(wf, atm_screen, **conf)
    
    # add static NCPAs
    if static_ncpa is True:
        filename = os.path.join(conf['input_dir'], conf['ncpa_screen_file'])
        ncpa_screen = fits.getdata(filename)
        ncpa_screen = np.nan_to_num(ncpa_screen)
        wf = static_ncpa(wf, ncpa_screen*1e-9, **conf)
    
    # add polishing errors
    if polish_error is True:
        filename = os.path.join(conf['input_dir'], conf['polish_screen_file'])
        polish_screen = fits.getdata(filename)
        polish_screen = np.nan_to_num(polish_screen)
        wf = static_ncpa(wf, polish_screen*1e-9, **conf)
    
    return wf
