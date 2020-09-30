from heeps.util.download_from_gdrive import extract_zip
import astropy.units as u
import os
import numpy as np
import proper
proper.print_it = False


def read_config(verbose=False, **update_conf):

    conf = dict(

    # =============================================================================
    #           Console and file management 
    # =============================================================================
    cpu_count = 1,                      # 1 = single core; None = max-1 cores
    send_to = None,                     # user's email, for notifications
    prefix = '',                        # for saved files: e.g. 'test_'
    headless = False,                   # true if running on a headless server
    gdrive_id = '1wj3onWQ9GVW-l8X58JMgAj-9TNqalKb-', # Google Drive ID
    # required directories for data (e.g. fits files)
    dir_current = os.getcwd(),
    dir_input = 'input_files',
    dir_output = 'output_files',
    dir_temp = 'temp_files',

    # =============================================================================
    #           Parameters for telescope
    # =============================================================================
    focal = 658.6,                      # focal distance in m
    pupil_img_size = 39.9988,           # pupil image in m (for PROPER)
    diam_ext = 36.905,                  # outer circular aperture in m
    diam_int = 11.213,                  # central obscuration in m
    spi_width = 0.5,                    # spider width in m
    spi_angles = [0, 60, 120],          # spider angles in deg
    file_pupil = 'pupils/ELT_pupil_1385.fits',# entrance pupil file
    # if no valid pupil file, pupil will be created with the following params:
    seg_width = 1.45,                   # segment width in m
    seg_gap = 0.004,                    # gap between segments in m
    seg_rms = 0,                        # rms of the reflectivity of all segments
    select_petal = None,
    npetals = 6,                        # number of petals
    # number of hexagonal segments per column (from left to right)
    seg_ny = np.array([10, 13, 16, 19, 22, 23, 24, 25, 26, 27, 28, 29, \
                      30, 31, 30, 31, 30, 31, 30, 31, 30, 31, 30, 29, \
                      28, 27, 26, 25, 24, 23, 22, 19, 16, 13, 10]),
    # coordinates of missing segments
    seg_missing = [],
#    seg_missing = [(-2,5),(-2,6),(-3,4),(-3,5),(-3,6),(-4,5),(-4,6)],


    # =============================================================================
    #           Parameters for observing modes and bands
    # =============================================================================
    # Types of various modes include (e.g. for METIS):
    #    1. mode = 'RAVC' for Ring Apodized Vortex Coronagraph
    #    2. mode = 'CVC'  for Classical Vortex Coronagraph
    #    3. mode = 'APP'  for Apodizing Phase Plate
    #    4. mode = 'CLC'  for Classical Lyot Coronagraph
    #    5. mode = 'ELT'  for no coronagraph (only telescope)
    # Default mode: L-band Ring Apodized Vortex
    mode = 'RAVC',                      # HCI mode
    band = 'L',                         # spectral band
    lam = 3.81e-6,                      # wavelength in m
    pscale = 5.47,                      # pixel scale in mas/pix
    ngrid = 1024,                       # number of pixels of the wavefront array
    npupil = 285,                       # number of pixels of the pupil
    ndet = 365,                         # size of the detector plane array
    hfov = 1,                           # (optional) half FOV in arcsec (updates ndet)
    mag = 5,                            # star magnitude at selected band
    mag_ref = 0,                        # reference magnitude for star and background fluxes
    flux_star = 8.999e+10,              # [e-/s] HCI-L long, mag 0 (Jan 21, 2020)
    flux_bckg = 8.878e+04,              # [e-/s/pix]
    ls_dRext = 0.03,                    # LS Rext undersize (% diam ext)
    ls_dRint = 0.05,                    # LS Rint oversize (% diam ext)
    ls_dRspi = 0.04,                    # LS spider oversize (% diam ext)
    ls_misalign = [0,0,0,0,0,0],        # Lyot stop misalignment
    vc_charge = 2,                      # vortex topological charge
    ravc_calc = True,                   # calculate RA params (Mawet2013)
    ravc_t = 0.76,                      # default RA trans (calc=False)
    ravc_r = 0.62,                      # default RA radius (calc=False)
    ravc_misalign = [0,0,0,0,0,0],      # RA misalignment
    clc_diam = 4,                       # CLC occulter diam in lam/D
    file_ravc_phase = '',               # ring apodizer files (optional)
    file_ravc_amp = '',
    file_app_phase = 'APP/app_phase_cut.fits',# apodizing phase plate files
    file_app_amp = '',
    app_strehl = 0.64,                   # APP Strehl ratio
    app_single_psf = 0.48,               # APP single PSF (4% leakage)
    student_distrib = True,              # use Student's distribution instead of Gaussian
    # Multiple spectral bands
    bands = ['L', 'M', 'N1', 'N2'],
    band_specs = {  
        'L': {'lam': 3.82e-6,
            'pscale': 5.47,
            'flux_star': 8.999e+10,                 # HCI-L long
            'flux_bckg': 8.878e+04,
            'ls_dRspi': 0.04},
        'M': {'lam': 4.8e-6,
            'pscale': 5.47,
            'flux_star': 2.452e+10,                 # CO ref
            'flux_bckg': 6.714e+05,
            'ls_dRspi': 0.04},
        'N1': {'lam': 8.65e-6,
            'pscale': 6.79,
            'flux_star': 3.684e+10,                 # GeoSnap N1
            'flux_bckg': 4.725e+07,
            'ls_dRspi': 0.05},
        'N2': {'lam': 11.25e-6, 
            'pscale': 6.79,
            'flux_star': 3.695e+10,                 # GeoSnap N2
            'flux_bckg': 1.122e+08,
            'ls_dRspi': 0.05},
        'N1a': {'lam': 8.65e-6, 
            'pscale': 10.78,
            'flux_star': 2.979e+10,                 # Aquarius N1
            'flux_bckg': 9.630e+07,
            'ls_dRspi': 0.05},
        'N2a': {'lam': 11.25e-6, 
            'pscale': 10.78,
            'flux_star': 2.823e+10,                 # Aquarius N2
            'flux_bckg': 2.142e+08,
            'ls_dRspi': 0.05}
        },
    # Multiple HCI modes
    modes = ['RAVC', 'CVC', 'CLC', 'APP', 'ELT'],
    mode_specs = {
        'RAVC': {'ls_dRint': 0.03},
        'CVC': {'ls_dRint': 0.05},
        'CLC': {'ls_dRint': 0.10}
        },

    # =============================================================================
    #           Parameters for wavefront
    # ============================================================================
    nframes = 20,                       # number of frames to crop the input data
    nstep = 1,                          # take 1 frame every nstep (cubesize = nframes/nstep)

    add_scao = False,                   # SCAO residuals cube
    file_scao = 'SCAO/cube_COMPASS_Oct2018_RandomWind_100screens.fits',
    
    add_ncpa = False,                   # NCPA phase screen
    file_ncpa = 'NCPA/NCPA_IMG_LMPP1-SCAO_PYR.fits', # IMG LM @ 3.7 um
    #file_ncpa_screen = 'NCPA_IMG_NQPP1-SCAO_DET.fits', # IMG NQ @ 10 um

    add_petal_piston = False,           # petal piston (island effect)
    file_petals = 'petals/petal%s_253.fits',# one petal
    rms_phase_sta = 35.9,               # static (nm)
    rms_phase_qlsf = 20,                # quasistatic low spatial freq (nm)
    rms_phase_qhsf = 20,                # quasistatic high spatial freq (nm)
    rms_phase_dyn = 40,                 # dynamic (nm)

    add_pointing_err = False,           # pointing errors
    rms_point_qsta = 0.4,               # quasistatic (mas)
    rms_point_dyn = 2,                  # dynamic (mas)

    add_apo_drift = False,              # apodizer drift
    ptv_drift = 0.01                    # (%)

    )                                   # end of default conf dict
 
    # =============================================================================
    #           Perform initialization actions
    # =============================================================================
    
    # update conf dictionary
    conf.update(**update_conf)    
    if verbose is True:
        print('Read config: band=%s, mode=%s'%(conf['band'], conf['mode']))
        print('\u203e'*12)
        print('   npupil=%s, pscale=%s mas, lam=%3.4E m'\
            %(conf['npupil'], conf['pscale'], conf['lam']))
        hfov = conf['ndet']/2*conf['pscale']/1e3
        hfov_lamD = hfov*u.arcsec.to('rad')/(conf['lam']/conf['diam_ext'])
        print('   ndet=%s, hfov=%s arcsec (%s lam/D)\n'%(conf['ndet'], \
            round(hfov,2), round(hfov_lamD,2)))
    
    # create directories
    conf['dir_current'] = os.path.normpath(os.path.expandvars(conf['dir_current']))
    conf['dir_input'] = os.path.join(conf['dir_current'], conf['dir_input'], '')
    conf['dir_output'] = os.path.join(conf['dir_current'], conf['dir_output'], '')
    conf['dir_temp'] = os.path.join(conf['dir_current'], conf['dir_temp'], '')
    os.makedirs(conf['dir_input'], exist_ok=True)
    os.makedirs(conf['dir_output'], exist_ok=True)
    os.makedirs(conf['dir_temp'], exist_ok=True)
    
    # create paths to fits files
    conf['file_pupil'] = os.path.join(conf['dir_input'], conf['file_pupil'])
    conf['file_scao'] = os.path.join(conf['dir_input'], conf['file_scao'])
    conf['file_ncpa'] = os.path.join(conf['dir_input'], conf['file_ncpa'])
    conf['file_petals'] = os.path.join(conf['dir_input'], conf['file_petals'])
    conf['file_ravc_phase'] = os.path.join(conf['dir_input'], conf['file_ravc_phase'])
    conf['file_ravc_amp'] = os.path.join(conf['dir_input'], conf['file_ravc_amp'])
    conf['file_app_phase'] = os.path.join(conf['dir_input'], conf['file_app_phase'])
    conf['file_app_amp'] = os.path.join(conf['dir_input'], conf['file_app_amp'])
    
    # downloading input files from Google Drive
    if not os.path.isfile(conf['file_scao']):
        print("Downloading input files from Google Drive to \n'%s'\n"%conf['dir_input'])
        extract_zip(conf['gdrive_id'], conf['dir_input'])
    
    # disable matplotlib display to run on a headless server
    if conf['headless'] is True:
        import matplotlib; matplotlib.use('agg')

    # sort alphabetically
    conf = {k: v for k, v in sorted(conf.items())}

    return conf
