import numpy as np
import proper
from astropy.io import fits

from island_effect_piston import island_effect_piston
from atmosphere import atmosphere
from static_ncpa import static_ncpa


def wavefront_abberations(wfo, conf,atm_screen, TILT):

    diam = conf['DIAM']
    gridsize = conf['GRIDSIZE']
    pixelsize = conf['PIXEL_SCALE'] 
    wavelength = conf['WAVELENGTH']
    lamda = proper.prop_get_wavelength(wfo)
    Debug = conf['DEBUG']
    Debug_print = conf['DEBUG_PRINT'] 

    Island_Piston = np.array(conf['ISLAND_PISTON'])
#    atm_screen = fits.getdata(input_dir + '/' + atm_screen)
    beam_ratio = pixelsize*4.85e-9/(wavelength/diam)
    npupil = np.ceil(gridsize*beam_ratio) # compute the pupil size --> has to be ODD (proper puts the center in the up right pixel next to the grid center)

    if npupil % 2 == 0:
        npupil = npupil +1

    if (isinstance(atm_screen, (list, tuple, np.ndarray)) == True) and (atm_screen.ndim >= 2): # when the atmosphere is present
        atmosphere(wfo, npupil, atm_screen, Debug_print, Debug)
    
    if (all(v == 0 for v in Island_Piston) == False): # when the piston is present
        island_effect_piston(wfo, npupil, Island_Piston, Debug_print, Debug)

    if conf['STATIC_NCPA'] == True:
        filename = conf['IMG_LM_SCAO']
        phase_screen = fits.getdata(conf['INPUT_DIR'] + filename)
        phase_screen = np.nan_to_num(phase_screen)
        phase_screen *= 10**-9
        static_ncpa(wfo, npupil, phase_screen)
    
    if (TILT.any != 0.): # when tip/tilt
        if (Debug_print==True):
            print("TILT: ", TILT)
            print("lamda: ", lamda)
        tiptilt=(np.multiply(TILT, lamda))/4 # translate the tip/tilt from lambda/D into RMS phase errors
        proper.prop_zernikes(wfo, [2,3], tiptilt) # 2-->xtilt, 3-->ytilt
    
    return wfo
