from heeps.util.download_from_gdrive import extract_zip
import collections
import proper
import os

conf = collections.OrderedDict()

# =============================================================================
#           Console and file management 
# =============================================================================
proper.print_it = False                     # allow intermediate prints from PROPER
conf['cpucount'] = 1                        # 1 = single core; None = max-1 cores
conf['send_to'] = None                      # email for sim end notification
conf['send_subject'] = 'HEEPS noreply'
conf['send_message'] = 'HEEPS simulation finished.'
conf['prefix'] = ''                         # for saved files: e.g. 'test_'

# required directories for data (e.g. fits files)
conf['current_dir'] = os.getcwd()           # default
#conf['current_dir'] = '$HOME/INSTRUMENTS/METIS/heeps_analysis'
#conf['current_dir'] = '/mnt/disk4tb/METIS/heeps_analysis'
conf['input_dir'] = 'input_files'
conf['output_dir'] = 'output_files'
conf['temp_dir'] = 'temp_files'

# create paths and directories
conf['current_dir'] = os.path.normpath(os.path.expandvars(conf['current_dir']))
conf['input_dir'] = os.path.join(conf['current_dir'], conf['input_dir'], '')
conf['output_dir'] = os.path.join(conf['current_dir'], conf['output_dir'], '')
conf['temp_dir'] = os.path.join(conf['current_dir'], conf['temp_dir'], '')
os.makedirs(conf['input_dir'], exist_ok=True)
os.makedirs(conf['output_dir'], exist_ok=True)
os.makedirs(conf['temp_dir'], exist_ok=True)


# =============================================================================
#           Define parameters for Telescope
# =============================================================================
conf['lam'] = 5e-6                         # wavelength in meters
conf['diam'] = 37                          # diameter of the telescope in meters
conf['R_obstr'] = 0.3                      # secondary obstruction in percentage
conf['spiders_width'] = 0.60               # width of spiders in meters
conf['spiders_angle'] = [0, 60, 120]       # angles of spiders
conf['N_mis_segments'] = 0                 # number of missing segments
conf['gridsize'] = 1024                    # (integer) grid size of the simulation array
conf['pscale'] = 5.21                      # pixel scale in mas/pix (e.g. METIS LM=5.21, NQ=10.78)
conf['focal'] = 658.6                      # focal distance in meters
# SCAO team currently uses circular pupil with secondary obstruction
conf['pupil_file'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits' # input pupil
# downloading input files from Google Drive
conf['gdriveID'] = '1wj3onWQ9GVW-l8X58JMgAj-9TNqalKb-'
if not os.path.isfile(os.path.join(conf['input_dir'], conf['pupil_file'])):
    print("Downloading input files from Google Drive to '%s'."%conf['input_dir'])
    extract_zip(conf['gdriveID'], conf['input_dir'])


# =============================================================================
#           Parameters for Wavefront aberrations
# =============================================================================
#conf['zernike'] = None          # zernike polynomials
conf['zern_inds'] = [2,3]       # zernike indices

#conf['petal_piston'] = None     # petal piston (island effect)
conf['npetals'] = 6             # number of petals in the pupil

#conf['atm_screen'] = None       # atmosphere SCAO phase screens
conf['atm_screen_file'] = 'cube_atm_100screens_Feb2018_RandomWind.fits'

#conf['ncpa_screen'] = None      # quasi-static NCPAs
conf['ncpa_screen_file'] = 'NCPA_IMG_LMPP1-SCAO_PYR.fits' # IMG LM @ 3.7 um
#conf['ncpa_screen_file'] = 'NCPA_IMG_NQPP1-SCAO_DET.fits' # IMG NQ @ 10 um


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
