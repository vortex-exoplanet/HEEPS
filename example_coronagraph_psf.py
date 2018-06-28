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


# =============================================================================
#                         ELT Pupil Plane
# =============================================================================

(npupil, wfo) = pupil(diam, gridsize, spiders_width, spiders_angle, pixelsize,  r_obstr, wavelength,
            pupil_file=pupil_file, missing_segments_number=0, Debug=Debug, Debug_print=Debug_print, prefix=prefix)  

# =============================================================================
#                      Wavefront abberations
# =============================================================================
TILT = np.array(((0,0),(.1,.1)))#np.array([0.,0.])

wavefront_abberations(wfo, npupil, atm_screen, NCPA,Island_Piston,TILT=TILT[1], Debug='False', 
          Debug_print='False', prefix='test')  


# =============================================================================
#       Coronagraph selection -- Classical Vortex / RAVC / APP --
# =============================================================================


#   (1) __________  Apodization (for RAVC and APP)  __________________________


RAVC = False

apodization(wfo, r_obstr, npupil, RAVC=RAVC, phase_apodizer_file=phase_apodizer_file, amplitude_apodizer_file=amplitude_apodizer_file, apodizer_misalignment=apodizer_misalignment, Debug_print=Debug_print)

#   (2)  __________  Coronagraph mask (for vortex or RAVC )  __________________


if (Vortex == True):
    vortex(wfo, charge, f_lens,diam, pixelsize, Debug_print = Debug_print)

##  (3)  _____________       lyot stop ________________________________________

lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, 
         LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug)


# =============================================================================
#       Detector plane
# =============================================================================

psf = detector(wfo,f_lens,nd)	


plt.figure()
plt.imshow(np.sqrt(psf))
plt.colorbar()

plt.show()

