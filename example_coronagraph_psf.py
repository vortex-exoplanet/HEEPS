#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a vortex PSF
# =============================================================================

# Required libraries
#import proper # library for propagation of wavefront
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits    # to read and write fits files

from simulation_config import *	# Loads default configuration for simulation (later will be changed to dictinary)
from heeps import * # loads all HEEPS scripts required for simuation


#  default simulation cofiguration defined in "simulation_config" can be overridden here
wavelength = 3.80*10**-6
charge = 2 # charge is modified here
TILT = np.array(((0,0),(.1,.1)))#np.array([0.,0.])


# =============================================================================
#                         ELT Pupil Plane
# =============================================================================

(npupil, wfo) = pupil(diam, gridsize, spiders_width, spiders_angle, pixelsize,  r_obstr, wavelength,
            pupil_file=pupil_file, missing_segments_number=0, Debug=Debug, Debug_print=Debug_print, prefix=prefix)  

# =============================================================================
#                      Wavefront abberations
# =============================================================================

wavefront_abberations(wfo, npupil, atm_screen, NCPA,Island_Piston,TILT=TILT[0], Debug='False', Debug_print='False', prefix='test')  


# =============================================================================
#       Coronagraph selection -- Vortex Classical (VC) / RAVC / APP --
# =============================================================================

coronagraphs(wfo, r_obstr,npupil, phase_apodizer_file,amplitude_apodizer_file,apodizer_misalignment,charge,
              f_lens,diam,LS_amplitude_apodizer_file,LS_misalignment,LS,LS_parameters,spiders_angle, LS_phase_apodizer_file, Debug_print,pixelsize, Debug,coronagraph_type= 'VC')


# =============================================================================
#       Detector plane
# =============================================================================

psf = detector(wfo,f_lens,nd)	


plt.figure()
plt.imshow(np.sqrt(psf))
plt.colorbar()

plt.show()

