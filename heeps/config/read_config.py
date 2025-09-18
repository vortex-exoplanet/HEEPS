from heeps.util.download_from_gdrive import extract_zip
import astropy.units as u
import os
import numpy as np
import proper
proper.print_it = False


def read_config(verbose=False, **update_conf):
    """Creates a configuration dictionary with all useful parameters
    set to their default values. The user can specify different values 
    either by editing the read_config file, or by passing a dictionary 
    as argument.

    Args:
        update_conf:
            user specified values

    Returns:
        conf (dict):
            configuration dictionary for end-to-end simulation 

    """

    conf = dict(
    
    # =============================================================================
    #           Console and file management
    # =============================================================================
    cpu_count = None,                   # 1 = single core; None = max cores
    send_to = None,                     # user's email, for notifications
    prefix = '',                        # for saved files: e.g. 'test_'
    headless = False,                   # true if running on a headless server
    gdrive_id = '1LKkDZMEVtAxF82QExsf8CsQvqqSRnRWa', # Google Drive ID September 2025 -- previously:'1eEpEM862b0u8L4_AtR9NGlLKSUI0euPb'
    # required directories for data (e.g. fits files)
    dir_current = '$HOME/heeps_metis',
    dir_input = 'input_files',
    dir_output = 'output_files',
    dir_temp = 'temp_files',

    # =============================================================================
    #           Parameters for telescope and optics
    # =============================================================================
    focal = 658.6,                      # focal distance in m
    pupil_img_size = 39.9988,           # pupil image in m (for PROPER)
    diam_nominal = 38.542,              # nominal diameter (for LS oversizing)
    diam_ext = 36.905,                  # effective outer circular aperture in m
    diam_int = 11.213,                  # effective central obscuration in m
    diam_vpm = 16.273,                  # vortex phase mask (19/45*diam_nominal)
    diam_lom = 16.273,                  # lyot occulting mask (19/45*diam_nominal)
    f_pupil = 'pupil/ELT_fullM1_feb2025.fits.tar.gz',   # entrance pupil file
    # if no valid pupil file, pupil will be created with the following params:
    spi_angles = [0,60,120,180,240,300], # regular spider angles in deg
    spi_width = 0.202,                   # regular spider width in m
    AP_angles = [180],                   # Asymetric Pupil spider angles in deg
    AP_width = 0.386,                   # AP spider width in m (AP_width > spi_width)
    AP_center = 1,                      # AP spider normalized center wrt pupil radius
    AP_length = 2,                      # AP spider normalized length wrt pupil radius
    add_seg = False,                    # if Tue, will add segments
    # number of hexagonal segments per column (from left to right)
    seg_ny = np.array([10, 13, 16, 19, 22, 23, 24, 25, 26, 27, 28, 29, \
                       30, 31, 30, 31, 30, 31, 30, 31, 30, 31, 30, 29, \
                       28, 27, 26, 25, 24, 23, 22, 19, 16, 13, 10]),
    seg_missing = [],                   # missing segments coord, e.g. [(-2,5),(-2,6)]
    seg_rint = 4.01,                    # number of black segments accross the central obstruction radius
    seg_width = 1.417,                  # segment width in m
    seg_gap = 0.004,                    # gap between segments in m
    seg_ptv = 0,                        # ptv uniform segment intensity loss
    select_petal = None,
    npetals = 6,                        # number of petals
    
    # =============================================================================
    #           Parameters for observing modes and bands
    # =============================================================================
    # Types of various modes include (e.g. for METIS):
    #    1. mode = 'RAVC' for Ring Apodized Vortex Coronagraph
    #    2. mode = 'CVC'  for Classical Vortex Coronagraph
    #    3. mode = 'APP'  for Apodizing Phase Plate
    #    4. mode = 'SPP'  for Shaped Pupil Plate
    #    5. mode = 'CLC'  for Classical Lyot Coronagraph
    #    6. mode = 'ELT'  for no coronagraph (only telescope)
    # Multiple HCI modes
    modes = ['RAVC', 'CVC', 'CLC', 'APP', 'SPP', 'ELT'],
    # Default mode: L-band Ring Apodized Vortex
    mode = 'RAVC',                      # HCI mode
    band = 'L',                         # spectral band
    lam = 3.81e-6,                      # wavelength in m
    pscale = 5.47,                      # pixel scale in mas/pix
    ngrid = 1024,                       # number of pixels of the wavefront array
    npupil = 285,                       # number of pixels of the pupil
    ndet = 365,                         # size of the detector plane array
    hfov = 1.1,                         # (optional) half FOV in arcsec (updates ndet)
    add_bckg = False,                   # true means background flux and photon noise are added 
    mag = 5,                            # star magnitude at selected band
    mag_ref = 0,                        # reference magnitude for star and background fluxes
    flux_star = 8.999e+10,              # [e-/s] HCI-L long, mag 0 (Jan 21, 2020)
    flux_bckg = 8.878e+04,              # [e-/s/pix]
    call_ScopeSim = False,              # true if interfacing ScopeSim
    duration = 3600,                    # duration of the ADI sequence in s
    dit = 0.3,                          # detector integration time in s
    lat = -24.59,                       # telescope latitude in deg (Armazones=-24.59 ,Paranal -24.63)
    dec = -5,                           # star declination in deg (e.g. 51 Eri -2.47)
    # Lyot stop: 1. auto_select_lyot based on mode and band. 2. lyot_stop_name. 3. f_lyot_stop. 4. Provide parameter to build Lyot stop
    auto_select_lyot = True,            # True if automatically select lyot stop from definition_lyotstop
    lyot_stop_name = '',                # leave empty if auto_select_lyot=True. Otherwise, see definition_lyotstops for list of available stops
    f_lyot_stop = '',                   # lyot stop file. Only if auto_select_lyot=False and lyot_stop_name is empty.
    ls_dRext = 0.0477,                  # LS Rext undersize (% diam ext)
    ls_dRint = 0.04,                    # LS Rint oversize (% diam ext)
    ls_dRspi = 0.0275,                  # LS spider oversize (% diam ext)
    force_sym = False,                  # Force LS spider symmetry ?
    ls_ext_circ = False,                # circular LS external diameter
    ls_int_circ = False,                # circular LS internal diameter
    ls_misalign = None,                 # constant lyot stop misalignment
    vc_charge = 2,                      # vortex topological charge
    vc_zoffset = 0,                     # vortex defocus in m (z axis)
    vc_chrom_leak = 2e-3,               # vortex chromatic leakage
    lt_dist_start = 0,                  # light trap downstream of VPM
    lt_diam_start = 0,                  # e.g. dist=188e-3, diam=80e-3
    lt_dist_end = 0,                    # light trap upstream of LS
    lt_diam_end = 0,                    # e.g. dist=36.5e-3, diam=50e-3
    add_cl_vort = False,                # add chromatic leakage at the vortex plane
    add_cl_det = False,                 # add chromatic leakage at the detector plane
    ravc_calc = False,                  # calculate RAP params (Mawet2013)
    ravc_t = 0.7909,                    # (calc=False) RAP amplitude trans
    ravc_r = 0.5190,                    # (calc=False) RAP radius wrt allglass
    apo_misalign = None,                # constant apodizer misalignment
    clc_diam = 80,                      # CLC occulter diam in mas
    f_oat = 'optics/vc/oat_L_RAVC.fits',   # vortex off-axis transmission
    f_vc_trans = 'optics/vc/agpm_trans.fits', # vortex transmittance
    f_app_trans = 'optics/app/metis_gvapp_tx.fits', # APP transmittance
    f_app_amp = 'optics/app/METIS_IMG_aper_rot.fits', # APP amplitude
    f_app_phase = 'optics/app/METIS_IMG_gvAPP_final_rot.fits', # APP phase
    f_spp_amp = 'optics/spp/SP_AP2_masked_rot_bin.fits', # Shaped Pupil Plate
    app_strehl = 0.6,                   # APP Strehl ratio
    app_single_psf = 0.48,              # APP single PSF (4% leakage)
    app_phase_ramp_params = {"offset": 0, "angle": 0},    # units [lambda/D] and [deg]
    student_distrib = True,             # use Student's distribution instead of Gaussian
    # Multiple spectral bands
    bands = ['L', 'M', 'N1', 'N2'],
    band_specs = {
        'L': {'lam': 3.81e-6,
            'pscale': 5.47,
            'flux_star': 8.999e+10,                 # HCI-L long
            'flux_bckg': 8.878e+04},
        'M': {'lam': 4.79e-6,
            'pscale': 5.47,
            'flux_star': 2.452e+10,                 # CO ref
            'flux_bckg': 6.714e+05},
        'N1': {'lam': 8.70e-6,
            'pscale': 6.79,
            'flux_star': 3.684e+10,                 # GeoSnap N1
            'flux_bckg': 4.725e+07},
        'N2': {'lam': 11.33e-6, 
            'pscale': 6.79,
            'flux_star': 3.695e+10,                 # GeoSnap N2
            'flux_bckg': 1.122e+08},
        'N1a': {'lam': 8.67e-6, 
            'pscale': 10.78,
            'flux_star': 2.979e+10,                 # Aquarius N1
            'flux_bckg': 9.630e+07},
        'N2a': {'lam': 11.21e-6, 
            'pscale': 10.78,
            'flux_star': 2.823e+10,                 # Aquarius N2
            'flux_bckg': 2.142e+08}
        },


    # =============================================================================
    #           Parameters for wavefront
    # ============================================================================
    nframes = 4,                      # number of frames to crop the input data
    nstep = 1,                          # take 1 frame every nstep (cubesize = nframes/nstep)
    nframes_avg = 2,                   # number of frames averaged for off-axis psf

    add_phase = True,                   # phase screens (SCAO residuals, NCPA, petal piston)
    f_phase = 'wavefront/COMPASS_201810_RandomWind_100screens_meters.fits',
    add_amp = False,                    # amplitude screens (Talbot effect)
    f_amp = 'wavefront/talbot_LM_rep_6_field_24.fits',
    ncpa_sta = 35.9,                    # static (nm rms)
    ncpa_qlsf = 20,                     # quasistatic low spatial freq (nm rms)
    ncpa_qhsf = 20,                     # quasistatic high spatial freq (nm rms)
    ncpa_dyn = 40,                      # dynamic (nm rms)

    add_point_err = False,              # pointing errors
    f_point_err = 'wavefront/point_all_3600s_300ms_L.fits',
    point_qsta = 0.4,                   # quasistatic (mas rms)
    point_dyn = 2,                      # dynamic (mas rms)
    zern = None,                        # constant Zernike, e.g. astigmatism 30 nm rms,
                                        # zern = [5, 30e-9]

    add_apo_drift = False,              # apodizer drift
    apo_drift = 0.02,                   # (% ptv)

    add_ls_drift = False,               # Lyot stop drift
    ls_drift = 0.02,                    # (% ptv)

    )                                   # end of default conf dict

    # =============================================================================
    #           Perform initialization actions
    # =============================================================================
    
    # update conf dictionary
    conf.update(**update_conf)
    if verbose is True:
        print('Default config: band=%s, mode=%s'%(conf['band'], conf['mode']))
        print('\u203e'*15)
        print('   npupil=%s, pscale=%s mas, lam=%3.4E m'\
            %(conf['npupil'], conf['pscale'], conf['lam']))
        hfov = conf['ndet']/2*conf['pscale']/1e3
        hfov_lamD = hfov*u.arcsec.to('rad')/(conf['lam']/conf['diam_ext'])
        print('   ndet=%s (-> hfov=%s arcsec, %s lam/D)\n'%(conf['ndet'],
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
    for filename in ['f_pupil', 'f_phase', 'f_amp', 'f_point_err', 'f_lyot_stop',
            'f_oat', 'f_vc_trans', 'f_app_trans', 'f_app_amp', 'f_app_phase', 'f_spp_amp']:
        conf[filename] = os.path.join(conf['dir_input'], conf[filename])
    
    # downloading input files from Google Drive
    verif_file = 'wavefront/COMPASS_201810_RandomWind_100screens_meters.fits'
    if not os.path.isfile(os.path.join(conf['dir_input'], verif_file)):
        print("Downloading input files from Google Drive to \n'%s'\n"%conf['dir_input'])
        extract_zip(conf['gdrive_id'], conf['dir_input'])
    
    # disable matplotlib display to run on a headless server
    if conf['headless'] is True:
        import matplotlib; matplotlib.use('agg')

    # sort alphabetically
    conf = {k: v for k, v in sorted(conf.items())}

    return conf