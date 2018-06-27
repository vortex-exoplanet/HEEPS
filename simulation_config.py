#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
import os

out_dir = str('./output_files/')
tmp_dir = str('./temp_files/')
input_dir = str('./input_files/')

os.makedirs(out_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)

# =============================================================================
#           Define parameters for Telescope
# =============================================================================
gridsize = 512                    # (integar) grid size of the simulation array
wavelength = 5.0*10**-6           # wavelength in micron
diam = 37.0                       # diameter of the telescope in meters 
r_obstr = 0.3                     # secondary obstruction in percentage??
spiders_width = 0.60              # width of ELT spiders in meters
spiders_angle = [0, 60.0, 120.0]  # angles of spiders 
missing_segments_number = 0       # number of mission segments

pupil_file = fits.getdata(input_dir + '/ELT_2048_37m_11m_5mas_nospiders_cut.fits')

prefix = 'test1'                # prefix for saving fits_files   
pixelsize = 5.0                 # plate scale in milli_arc_sec/pixel
Debug=True                      # various fits files for coronagraphic propagation is saved in out_dir 
Debug_print = False             # prints various values   
f_lens = 658.6                  # float, meters, focal distance


# =============================================================================
#           Parameters for Wavefront abberations
# =============================================================================
TILT = np.array([0.,0.])
Island_Piston = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
NCPA = 0
atm_screen = fits.getdata(input_dir+'metis_370P_35L_HCI_Feb18_rwf8160_cut.fits')
#atm_screen = 0
#atm_screen = fits.getdata(input_dir+'cube_atm_1000screens_Feb2018_RandomWind.fits')[10]



# =============================================================================
#           Parameters for Coronagraphs
# =============================================================================
""" Types of coronagraphs include:
    1. Classical vortex: 
        to use this coronagraph; vortex = True and define charge 
    2. Ring apodized vortex coronagraph (RAVC):
        to use RAVC = True with the parameters for classical vortex 
    3. Apodizing phase plate (APP):
        to use APP input phase_apodizer_file as fits and 
        vortex = False, RAVC = False
"""

Vortex = True
charge = 2
RAVC = False
APP = False

if (APP==True):
    Vortex = False
    RAVC = False
    phase_apodizer_file = fits.getdata(input_dir+'app_phase_cut.fits')
else:
    phase_apodizer_file = 0


amplitude_apodizer_file = 0 
apodizer_misalignment = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

LS = True
LS_parameters = [0.98, 0.03, 1.1]
LS_phase_apodizer_file = 0
LS_amplitude_apodizer_file = 0
LS_misalignment = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


# =============================================================================
#           Detector plane
# =============================================================================

nd = 400                # final image size





