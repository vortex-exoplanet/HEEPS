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
n = gridsize

# =============================================================================
#                         ELT Pupil Plane
# =============================================================================


# =============================================================================
#                      Wavefront abberations
# =============================================================================

(npupil, wfo1) = pupil(diam, gridsize, spiders_width, spiders_angle, pixelsize,  r_obstr, wavelength,
        pupil_file=pupil_file, missing_segments_number=0, Debug=Debug, Debug_print=Debug_print, prefix=prefix) 

if (atm_screen.ndim == 3) or (TILT.ndim == 2) or (LS_misalignment.ndim == 2) or (apodizer_misalignment.ndim == 2) or (Island_Piston.ndim == 2) or (NCPA.ndim == 3):
    print('Cube')
    
    if (atm_screen.ndim == 3):
        length_cube = atm_screen.shape[0]
    if (TILT.ndim == 2):
        length_cube = TILT.shape[0]
    if (LS_misalignment.ndim == 2):
        length_cube = LS_misalignment.shape[0]
    if (apodizer_misalignment.ndim == 2):
        length_cube = apodizer_misalignment.shape[0]
    if (Island_Piston.ndim == 2):
        length_cube = Island_Piston.shape[0]
    if (NCPA.ndim == 3):
        length_cube = NCPA.shape[0]
    
    wfo = deepcopy(wfo1)
    psf_Coro = np.zeros((length_cube,nd,nd))

    
    for iter in range(0, length_cube):
        print('iter: ', iter)

        if ((isinstance(atm_screen, (list, tuple, np.ndarray)) == True)):
            if (atm_screen.ndim == 3):
                atm_screen_iter = atm_screen[iter,:,:]
            else:
                atm_screen_iter = atm_screen
        if (TILT.ndim == 2):
            TILT_iter = TILT[iter,:]
            print('TILT: ', TILT_iter)
        else:
            TILT_iter = TILT
        if (LS_misalignment.ndim == 2):
            LS_misalignment_iter =  LS_misalignment[iter,:]
        else:
            LS_misalignment_iter =  LS_misalignment
        if (apodizer_misalignment.ndim == 2):
            apodizer_misalignment_iter = apodizer_misalignment[iter,:]
        else:
            apodizer_misalignment_iter = apodizer_misalignment
        if (Island_Piston.ndim == 2):
            Island_Piston_iter = Island_Piston[iter,:]
        else:
            Island_Piston_iter = Island_Piston
        if (isinstance(NCPA, (list, tuple, np.ndarray)) == True):
            if (NCPA.ndim == 3):
                NCPA_iter = NCPA[iter,:,:]
            else:
                NCPA_iter = NCPA
             

        wavefront_abberations(wfo, npupil, atm_screen_iter, NCPA,Island_Piston,TILT=TILT_iter, Debug='False',Debug_print='False', prefix='test')  
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
        psf_Coro[iter,:,:] = psf

        
    fits.writeto(prefix+'_psf_cube_Coro_nonorm.fits', psf_Coro, overwrite=True)
        
plt.figure()
plt.imshow(np.sqrt(psf_Coro[1]))
plt.colorbar()

plt.show()

