#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
import os, sys
import proper
# to create directories for inputing and storing fits files
out_dir = str('./output_files/')
tmp_dir = str('./temp_files/')
input_dir = str('./input_files/')

os.makedirs(out_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)



# =============================================================================
#           Define parameters for Telescope
# =============================================================================
gridsize = 1024                    # (integar) grid size of the simulation array
wavelength = 5.0*10**-6           # wavelength in micron
diam = 37.0                       # diameter of the telescope in meters 
r_obstr = 0.3                     # secondary obstruction in percentage??
spiders_width = 0.60              # width of ELT spiders in meters
spiders_angle = [0, 60.0, 120.0]  # angles of spiders 
missing_segments_number = 0       # number of mission segments

# to input pupil as a fits file, 
# SCAO team currently uses circular pupil with secondary obstruction
pupil_file = fits.getdata(input_dir + '/ELT_2048_37m_11m_5mas_nospiders_cut.fits')

prefix = 'test'                 # prefix for saving fits_files   
pixelsize = 5.0                 # plate scale in milli_arc_sec/pixel
f_lens = 658.6                  # float, meters, focal distance

Debug = False                   # various fits files for coronagraphic propagation is saved in out_dir 
Debug_print = False             # prints various values   


# =============================================================================
#           Parameters for Wavefront abberations
# =============================================================================
TILT = np.array([0.0,0.]) 
#TILT = np.random.randn(10,2) # creates an array of (n,2) tip/tilt values

atm_screen = np.array([0.0])        # No phase screen 
#atm_screen = fits.getdata(input_dir+'metis_370P_35L_HCI_Feb18_rwf8160_cut.fits') # Single phase screen

atm_screen_cube = fits.getdata(input_dir+'cube_atm_1000screens_Feb2018_RandomWind.fits') # multi-cube phase screen


Island_Piston = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Not tested yet
NCPA = np.array([0.])   # Not tested yet


# =============================================================================
#           Parameters for Coronagraphs
# =============================================================================
""" Types of coronagraphs include:
    1. Vortex coronagraph (VC)
    2. Ring apodized vortex coronagraph (RAVC):
    3. Apodizing phase plate (APP):
(a) By changing the "coronagraph_type" to "VC/RAVC/APP" coronagraphs can be selcted. 
(b) If the "coronagraph_type" is "OFFAXIS" a non-coronagraphic PSF with lyot-stop is generated
(c) If the "coronagraph_type" is anything except above keywords a normal ELT PSF is generated

"""

coronagraph_type = 'VC'
charge = 2 # For vortex coronagraph

# Apodizer phase file for APP coronagraph
phase_apodizer_file = fits.getdata(input_dir+'app_phase_cut.fits')
apodizer_misalignment = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

amplitude_apodizer_file = 0  

LS = True
LS_parameters = [0.98, 0.03, 1.1]
LS_misalignment = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])



LS_phase_apodizer_file = 0      # Not tested yet
LS_amplitude_apodizer_file = 0  # Not tested yet

# =============================================================================
#           Detector plane
# =============================================================================

nd = 512                # final image size


proper.print_it = False     # To suppress the printing of intermediate steps by PROPER routine



