import numpy as np
import proper

from island_effect_piston import island_effect_piston
from atmosphere import atmosphere

def wavefront_abberations(wfo, npupil, atm_screen,NCPA,Island_Piston, TILT=[0,0], Debug='False', 
          Debug_print='False', prefix='test'):
    
    lamda = proper.prop_get_wavelength(wfo)
    
    if (isinstance(atm_screen, (list, tuple, np.ndarray)) == True) and (atm_screen.ndim >= 2): # when the atmosphere is present
        atmosphere(wfo, npupil, atm_screen, Debug_print, Debug)
    
    if (all(v == 0 for v in Island_Piston) == False): # when the piston is present
        island_effect_piston(wfo, npupil, Island_Piston, Debug_print, Debug)
    
    if (TILT.any != 0.): # when tip/tilt
        if (Debug_print==True):
            print("TILT: ", TILT)
            print("lamda: ", lamda)
        tiptilt=(np.multiply(TILT, lamda))/4 # translate the tip/tilt from lambda/D into RMS phase errors
        proper.prop_zernikes(wfo, [2,3], tiptilt) # 2-->xtilt, 3-->ytilt
    
    return wfo
