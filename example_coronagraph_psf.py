#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic PSF
# =============================================================================

# Required libraries
import proper # library for propagation of wavefront
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits    # to read and write fits files

from simulation_config import *	# Loads default configuration for simulation 
from heeps import * # loads all HEEPS scripts required for simuation

"""
Default simulation cofiguration defined in "simulation_config.py" can be overridden here
"""
wavelength = 3.80*10**-6
charge = 2 # charge is modified here


# =============================================================================
#                         ELT Pupil Plane
# =============================================================================

(npupil, wfo) = pupil(diam, gridsize, spiders_width, spiders_angle, pixelsize, 
    r_obstr, wavelength, pupil_file=pupil_file, missing_segments_number=0, 
    Debug=True, Debug_print=Debug_print, prefix=prefix)  

# =============================================================================
#                      Wavefront abberations
# =============================================================================

wavefront_abberations(wfo, npupil, atm_screen, NCPA,Island_Piston,TILT=TILT, 
                      Debug='False', Debug_print='False', prefix='test')  

# =============================================================================
#       Coronagraph selection -- Vortex Classical (VC) / RAVC / APP --
# =============================================================================
"""
1. By changing the "coronagraph_type" to "VC/RAVC/APP" coronagraphs can be selcted. 
2. If the input is "None" a non-coronagraphic PSF with lyot-stop is generated
3. If the input is anything except above keywords a normal PSF is generated
"""

coronagraph_type = 'APP'

coronagraphs(wfo, r_obstr,npupil, phase_apodizer_file,amplitude_apodizer_file,
    apodizer_misalignment,charge,f_lens,diam,LS_amplitude_apodizer_file,LS_misalignment,
    LS,LS_parameters,spiders_angle, LS_phase_apodizer_file, Debug_print,pixelsize, 
    Debug,coronagraph_type= coronagraph_type)

# =============================================================================
#       Detector plane
# =============================================================================
psf = detector(wfo,f_lens,nd,coronagraph_type,prefix,Debug=True)	


plt.figure()
plt.imshow(np.sqrt(psf))
plt.colorbar()

plt.show()

