import numpy as np
import proper
from astropy.io import fits

from .island_effect_piston import island_effect_piston
from .atmosphere import atmosphere
from .static_ncpa import static_ncpa


def wavefront_abberations(wfo, AO_residuals=None, Island_Piston=None, tip_tilt=None,
            STATIC_NCPA=False, **conf):
    
    # add AO residuals
    if AO_residuals is not None:
        atmosphere(wfo, AO_residuals, **conf)
    
    # add island effect
    if Island_Piston is not None:
        island_effect_piston(wfo, npupil, Island_Piston, Debug_print, Debug)
    
    # add tip/tilt
    if tip_tilt is not None:
        lam = conf['WAVELENGTH'] # in meters
        tip_tilt = np.array(tip_tilt)*lam/4 # translate the tip/tilt from lambda/D into RMS phase errors
        proper.prop_zernikes(wfo, [2,3], tip_tilt) # 2-->xtilt, 3-->ytip?
    
    # add static NCPAs
    if STATIC_NCPA == True:
        filename = conf['IMG_LM_SCAO']
        phase_screen = fits.getdata(conf['INPUT_DIR'] + filename)
        phase_screen = np.nan_to_num(phase_screen)
        phase_screen *= 10**-9          # scale the wavelenth to nm
        static_ncpa(wfo, npupil, phase_screen)
    
    return wfo
