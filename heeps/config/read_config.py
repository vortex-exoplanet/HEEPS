import collections
import os
import proper
import numpy as np

conf = collections.OrderedDict()

# folder containing INPUT, OUTUT, and TEMP directories for data (e.g. fits files)
conf['data_folder'] = '.'

# Creates required directories for simulation 
conf['OUT_DIR'] = os.path.join(conf['data_folder'], 'output_files', '')
conf['TMP_DIR'] = os.path.join(conf['data_folder'], 'temp_files', '')
conf['INPUT_DIR'] = os.path.join(conf['data_folder'], 'input_files', '')
os.makedirs(conf['OUT_DIR'] , exist_ok=True)
os.makedirs(conf['TMP_DIR'], exist_ok=True)
os.makedirs(conf['INPUT_DIR'], exist_ok=True)


# =============================================================================
#           Define parameters for Telescope
# =============================================================================
conf['GRIDSIZE'] = 1024                    # (integar) grid size of the simulation array
conf['WAVELENGTH'] = 5.0*10**-6            # wavelength in meters
conf['DIAM'] = 37.0                        # diameter of the telescope in meters 
conf['R_OBSTR'] = 0.3                      # secondary obstruction in percentage??
conf['SPIDERS_WIDTH'] = 0.60               # width of ELT spiders in meters
conf['SPIDERS_ANGLE'] = [0, 60.0, 120.0]   # angles of spiders 
conf['MIS_SEGMENTS_NU'] = 0                # number of missing segments

# to input pupil as a fits file, 
# SCAO team currently uses circular pupil with secondary obstruction
conf['PUPIL_FILE'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
conf['PREFIX'] = 'test'                 # prefix for saving fits_files   
conf['PIXEL_SCALE'] = 5.0               # plate scale in milli_arc_sec/pixel
conf['F_LENS'] = 658.6                  # float, meters, focal distance
conf['DEBUG'] = False                   # various fits files for coronagraphic propagation is saved in out_dir 
conf['DEBUG_PRINT'] = False             # prints various values   

# =============================================================================
#           Parameters for Wavefront abberations
# =============================================================================
conf['TILT_2D'] = [0.0, 0.0]
conf['TILT_CUBE'] = 10                  # creates an array of (n,2) tip/tilt values
conf['ATM_SCREEN_NO'] = 0               # No phase screen 
conf['ATM_SCREEN_2D'] = 'metis_370P_35L_HCI_Feb18_rwf8160_cut.fits'   # Single phase screen
conf['ATM_SCREEN_CUBE'] = 'cube_atm_1000screens_Feb2018_RandomWind.fits' # 1000 phase screens
#conf['ATM_SCREEN_CUBE'] = 'atm_cube_100ms.fits' # 6000 phase screens
conf['GDRIVE_ID'] = '1AUtELRfn_xjnbsMM_SJG6W0c26zyzqNH'   # google drive id linked to the fits file 

conf['ISLAND_PISTON'] = None
conf['STATIC_NCPA'] = True
conf['IMG_LM_SCAO'] = 'NCPA_IMG_LMPP1-SCAO_DET.fits'          # NCPA phase screen @ 3.7 um
conf['IMG_NQ_SCAO'] = 'NCPA_IMG_NQPP1-SCAO_DET.fits'          # NCPA phase screen @ 10 um



# =============================================================================
#           Parameters for Coronagraphs
# =============================================================================
""" Types of various modes include:
    1. Vortex coronagraph (VC), to use put conf['MODE'] = 'VC'
    2. Ring apodized vortex coronagraph (RAVC), to use put conf['MODE'] = 'RAVC'
    3. Apodizing phase plate (APP), to use put conf['MODE'] = 'APP'
    4. No coronagraph, just Lyot-stop, to use put conf['MODE'] = 'OFFAXIS'    
    5. No coronagraph, but Ring apodizer and LS present, to use put conf['MODE'] = 'MASK'
    6. If conf['MODE'] = anything except above keywords ELT psf is generated    
"""

conf['MODE'] = 'VC'                     # default is vortex coronagraph   
conf['CHARGE'] = 2                      # default is charge 2 (AGPM)

conf['PHASE_APODIZER_FILE'] = 'app_phase_cut.fits'
conf['APODIZER_MIS_ALIGN'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
conf['AMP_APODIZER'] = 0

conf['LYOT_STOP'] = True
conf['LS_PARA'] = [0.98, 0.03, 1.1]
conf['LS_MIS_ALIGN'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
conf['DIAM_CL'] = 4     # classical lyot diam in lam/D

# =============================================================================
#           Detector plane
# =============================================================================
conf['N_D'] = 512                # final image size


collections.OrderedDict(sorted(conf.items()))



proper.print_it = False   # To suppress the printing of intermediate steps by PROPER routine
