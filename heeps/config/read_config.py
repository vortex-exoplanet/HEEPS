import collections
import os
import proper
import zipfile
from heeps.config import download_from_gdrive

conf = collections.OrderedDict()

# =============================================================================
#           Console and file management 
# =============================================================================

# allow PROPER to print intermediate steps
proper.print_it = False

# optional prefix for saved files: e.g. 'test_'
conf['prefix'] = ''

# create required directories for data (e.g. fits files)
conf['current_dir'] = '.'
conf['output_dir'] = os.path.join(conf['current_dir'], 'output_files', '')
conf['temp_dir'] = os.path.join(conf['current_dir'], 'temp_files', '')
conf['input_dir'] = os.path.join(conf['current_dir'], 'input_files', '')
os.makedirs(conf['output_dir'], exist_ok=True)
os.makedirs(conf['temp_dir'], exist_ok=True)
os.makedirs(conf['input_dir'], exist_ok=True)

# Google Drive parameters
conf['testfile'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
conf['gdriveID'] = '1YR4G_8E7TpznTxumQA1zD4z_v4V5ZlF2'
conf['gdriveZip'] = 'fits_files.zip'
# if test file is missing, start downloading from Google drive
if not os.path.isfile(os.path.join(conf['input_dir'], conf['testfile'])): 
    download_from_gdrive(conf['gdriveID'], conf['current_dir'], conf['gdriveZip'])
    with zipfile.ZipFile(conf['gdriveZip'],'r') as zip_ref:
        zip_ref.extractall(conf['current_dir'])
    os.remove(conf['gdriveZip'])

# server parameters, multiprocessing
conf['cpucount'] = 1                       # 1 = single core; None = max-1 cores
conf['send_subject'] = 'HEEPS noreply'
conf['send_message'] = 'HEEPS simulation finished.'
conf['send_to'] = None

# =============================================================================
#           Define parameters for Telescope
# =============================================================================
conf['lam'] = 5e-6                         # wavelength in meters
conf['diam'] = 37.0                        # diameter of the telescope in meters
conf['R_obstr'] = 0.3                      # secondary obstruction in percentage
conf['spiders_width'] = 0.60               # width of spiders in meters
conf['spiders_angle'] = [0, 60.0, 120.0]   # angles of spiders
conf['N_mis_segments'] = 0                 # number of missing segments
conf['gridsize'] = 1024                    # (integer) grid size of the simulation array
conf['pscale'] = 5.21                      # pixel scale in mas/pix (e.g. METIS LM=5.21, NQ=10.78)
conf['focal'] = 658.6                      # focal distance in meters
# SCAO team currently uses circular pupil with secondary obstruction
conf['pupil_file'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits' # input pupil


# =============================================================================
#           Parameters for Wavefront aberrations
# =============================================================================
conf['tip_tilt'] = [0., 0.]
conf['atm_screen_cube'] = 'cube_atm_1000screens_Feb2018_RandomWind.fits'
# google drive ID linked to the SCAO cube fits file 
conf['atm_screen_gdriveID'] = '1AUtELRfn_xjnbsMM_SJG6W0c26zyzqNH'

conf['petal_piston'] = None
conf['static_ncpa'] = False
conf['polish_error'] = False
conf['ncpa_screen'] = 'NCPA_IMG_LMPP1-SCAO_PYR.fits' # IMG LM @ 3.7 um
#conf['ncpa_screen'] = 'NCPA_IMG_NQPP1-SCAO_DET.fits' # IMG NQ @ 10 um


# =============================================================================
#           Parameters for Coronagraphs
# =============================================================================
""" Types of various modes include (e.g. for METIS):
    0. conf['mode'] = 'ELT'  for no coronagraph (only telescope)
    1. conf['mode'] = 'RAVC' for Ring Apodized Vortex Coronagraph
    2. conf['mode'] = 'CVC'  for Classical Vortex Coronagraph
    3. conf['mode'] = 'APP'  for Apodizing Phase Plate
    4. conf['mode'] = 'CLC'  for Classical Lyot Coronagraph
An off-axis PSF can be obtained by switching conf['onaxis'] to False,
thereby decentering the focal plane mask (if any).
"""
conf['onaxis'] = True                           # True for mask centered
conf['mode'] = 'CVC'                            # default is vortex coronagraph
conf['VC_charge'] = 2                           # default is charge 2 (AGPM)
conf['RAVC_misalign'] = [0,0,0,0,0,0]           # ring apodizer misalignment
conf['LS_misalign'] = [0,0,0,0,0,0]             # Lyot stop misalignment
conf['LS_params'] = [0.98, 0.03, 1.1]           # R_out(%), dR_in(%), LS spider width(m)
conf['APP_phase_file'] = 'app_phase_cut.fits'   # apodizing phase plate files
conf['APP_amp_file'] = ''
conf['CLC_diam'] = 4                            # CLC occulter diam in lam/D


# =============================================================================
#           Detector plane
# =============================================================================
conf['ndet'] = 512                              # final image size on detector

collections.OrderedDict(sorted(conf.items()))
