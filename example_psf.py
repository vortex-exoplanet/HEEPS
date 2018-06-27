#!/usr/bin/env python3

# =============================================================================
#       Example file for creating simple PSF
# =============================================================================

# Required libraries
import proper # library for propagation of wavefront
import matplotlib.pyplot as plt # for plotting simulated PSFs
import os

from simulation_config import *	# Loads default configuration for simulation (later will be changed to dictinary)

from heeps import pupil, detector, wavefront_abberations # loads HEEPS scripts required for simuation


out_dir = str('./output_files/')
tmp_dir = str('./temp_files/')
input_dir = str('./input_files/')

os.makedirs(out_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)



wavelength = 5*10**-6


# =============================================================================
#       ELT Pupil Plane
# =============================================================================

wfo = pupil(diam, gridsize, spiders_width, spiders_angle, pixelsize,  r_obstr, wavelength,
      missing_segments_number=0,  ELT_circ='True', Debug=Debug, Debug_print=Debug_print, prefix=prefix)

# =============================================================================
#       Wavefront abberations
# =============================================================================

TILT = np.array([0,10])

wavefront_abberations(wfo, gridsize, atm_screen, TILT, NCPA,Island_Piston, Debug='False', Debug_print='False', prefix='test')   


# =============================================================================
#       Detector plane
# =============================================================================
nd = 512

n = gridsize

psf = detector(wfo,f_lens)	

psf = psf[int(n/2-nd/2):int(n/2+nd/2),int(n/2-nd/2):int(n/2+nd/2)]


plt.figure()
plt.imshow((psf)**0.25)
plt.colorbar()



plt.show()

