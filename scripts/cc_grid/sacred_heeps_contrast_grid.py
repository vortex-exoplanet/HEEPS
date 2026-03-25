from sacred import Experiment, Ingredient
from sacred.observers import FileStorageObserver
import sys
import time
import os
import heeps
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments
import matplotlib.pyplot as plt
import astropy.io.fits as fits
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from create_phase_map_all import PhaseCubeGenerator
from fingerprint import writeToFile, getFingerprint

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import threadpoolctl
threadpoolctl.threadpool_limits(limits=1, user_api='blas')
print(threadpoolctl.threadpool_info())

import multiprocessing as mpro
try:
    mpro.set_start_method('fork')
except RuntimeError:
    print('[WARN] context has already been set')

HOMEDIR='/home/ulg/PSILab/gorban/'

# Initialize the experiment
# ex = Experiment("contrast_curve_simulation")
ex_name = "2026_contrast_curve_tests"
ex = Experiment("2026_contrast_curve_tests")

# ex.observers.append(FileStorageObserver.create("runs"))  # Save runs to a local directory
# db_path = os.path.dirname(__file__) + '/sacred_db/'

# db_path='/mnt/disk20tb/METIS/DATABASE_SACRED/' + ex_name + '/'
db_path='/globalscratch/ulg/PSILab/gorban/DATABASE_SACRED/' + ex_name + '/'
ex.observers.append(FileStorageObserver(db_path,
                                        priority=31))

# --- Define Ingredients ---
ncpa_ingredient = Ingredient("ncpa")

# Add ingredients to the experiment
ex.ingredients = [
    ncpa_ingredient,
]

# --- Default Configurations for Ingredients ---

@ncpa_ingredient.config
def default_ncpa():
    frequency = 10  # Hz (default)
    lag = 2
    gain_I=0.33
    nmodes=20

# ---------------------------------------------------------------------------
# Module-level constants — immutable, not logged by Sacred
# ---------------------------------------------------------------------------
NPUPIL = {'L': 285, 'M': 227, 'N1': 155, 'N2': 119}
HFOV = {'L': 0.8, 'M': 0.8, 'N1': 1.3, 'N2': 1.3}

PHASE_DIR  = 'wavefront/dfull/2026_grid/'
WV_DIR     = 'wavefront/wv/2026_grid/'
OAT_DIR    = 'optics/vc/'

_ncpa_dir   = 'wavefront/cbw/20260127/ncpa/'
CBW_FILENAMES = {
    'L':  _ncpa_dir + 'L_rep_1.fits',
    'M':  _ncpa_dir + 'M_rep_1.fits',
    'N1': _ncpa_dir + 'N1_rep_1.fits',
    'N2': _ncpa_dir + 'N2_rep_1.fits',
}

_talbot_dir = 'wavefront/cbw/20260127/talbot/'
TALBOT_FILENAMES = {
    'L':  _talbot_dir + 'L_rep_1.fits',
    'M':  _talbot_dir + 'M_rep_1.fits',
    'N1': _talbot_dir + 'N1_rep_1.fits',
    'N2': _talbot_dir + 'N2_rep_1.fits',
}

PUP_FILENAMES = {
    band: PHASE_DIR + f'mask4heeps_Telescope_Pupil_{npupil}.fits'
    for band, npupil in NPUPIL.items()
}

# RMS water-vapour OPD [nm]  (E-TNT-ULG-MET-0044 v0-1, Table 3-1)
WV_TABLE = {
    'L':  {'Q1': 14,  'Q2': 21,  'Q3': 29,  'Q4':  40, 'median':  24},
    'M':  {'Q1': 42,  'Q2': 60,  'Q3': 83,  'Q4': 115, 'median':  71},
    'N1': {'Q1': 79,  'Q2': 114, 'Q3': 158, 'Q4': 220, 'median': 136},
    'N2': {'Q1': 178, 'Q2': 257, 'Q3': 356, 'Q4': 495, 'median': 305},
}

# Magnitude threshold below which 10Hz NCPA correction is applied
# for greater or equal values, uses NCPA_FREQ_FALLBACK
NCPA_MAG_THRES = {
    'L':  {'CVC': 6.5, 'RAVC': 4.5},
    'M':  {'CVC': 3.5},
    'N1': {'CVC': 2.5},
    'N2': {'CVC': 1.5},
}
# NCPA loop frequency used when the star is too faint for the 10Hz nominal cadence
NCPA_FREQ_FALLBACK = 1  # Hz

# SCAO K-magnitude step-up thresholds: list of (mag_threshold, K) in descending order.
# First threshold that magnitude meets or exceeds wins; falls back to default K=6.
SCAO_K_THRES = [(8.5, 9), (7.5, 8)]


# QACITS residual tip-tilt error in mas
# valid at all magnitude because ALF sets the transition from 10 to 1Hz (not QACITS)
QACITS_ERROR_MAS = {'L':0.5, 'M': 0.5, 'N1':1, 'N2':1.5}
D_TEL = 36.905 # used for conversion from mas to nm rms

# MASK TAG used to select the file and get the ALF sensor error
MASK_TAGS={'L':  {'CVC': 'CLS-LM', 'RAVC': 'RLS-LM'},
           'M':  {'CVC': 'CLS-LM'},
           'N1': {'CVC': 'CLS-N'},
           'N2': {'CVC': 'CLS-N'}}
# ---------------------------------------------------------------------------


# --- Main Experiment Config ---
@ex.config
def default_config():
    tag=None        # this is used as prefix for HEEPS files (onaxis offaxis)
    mode = "CVC"  # "CVC" or "RAVC"
    band = "L"  # Default band
    nframes=36000
    duration = 3600
    dit = 0.1
    # Default magnitude -- get mapped to heeps config 'mag' in ``run_simulation``
    magnitude = 6.0   
    
    dry_run = False             # if True, skip all heavy computation (phase cube, propagation, contrast curves)
    do_f_phase = True           # whether to (re-)generate the combined phase cube
    do_propagation = True       # whether to run HEEPS wavefront propagation
    do_contrast_curves = True   # whether to compute contrast curves

    dir_current = HOMEDIR+'heeps_metis/' #'/home/gorban/heeps_metis/'
    dir_input = os.path.join(dir_current, 'input_files')

    cpu_count=10
    
    # Atmosphere parameters
    seeing = "Q2"  # Percentile Q1, Q2, Q3, Q4 = [12.5%, 37.5%, 62.5%, 87.5%]
    
    # SCAO configuration
    scao_K = 6  # Default SCAO mode (overridden based on magnitude in derived_config)

    # Other HEEPS specific variables
    nstep=1
    nframes_avg=100
    dec=-5
    add_phase=True
    add_amp=True

    select_lyot='auto'
    f_lyot_stop=''

    # normally set to False in HEEPS config by setting them explicitly here:
    add_cl_vort=False
    add_point_err=False
    add_seg=False
    add_apo_drift=False


# --- Derived Configuration (logged per-run by Sacred) ---
@ex.config
def derived_config(band, magnitude, mode, seeing, ncpa, do_f_phase,
                   duration, dit, dir_input, dir_current, scao_K):
    """
    Compute all derived scalars and input file paths from the base config and
    ingredient configs.  Everything assigned here is automatically recorded in
    Sacred's config.json for each run, making the grid fully reproducible.
    """

    # SCAO K-magnitude: override default (6) if magnitude crosses thresholds
    for _thresh, _k in SCAO_K_THRES:
        if magnitude >= _thresh:
            scao_K = _k
            break

    # RMS water-vapour OPD [nm]
    wv_rms = WV_TABLE[band][seeing]

    # NCPA AO-loop frequency [Hz]
    ncpa_freq = ncpa['frequency'] if magnitude <= NCPA_MAG_THRES[band][mode] else NCPA_FREQ_FALLBACK

    # Input file paths
    npupil   = NPUPIL[band]
    hfov     = HFOV[band]
    f_pupil  = PUP_FILENAMES[band]
    f_cbw    = CBW_FILENAMES[band]
    f_amp = TALBOT_FILENAMES[band]
    f_oat    = OAT_DIR    + f'oat_{band}_{mode}.fits'
    f_wv     = (WV_DIR     + f'cube_WV_20260225_3600_100ms_'
                          + f'Kmag2_0piston_meters_scao_only_{npupil}.fits')
    f_scao   = (PHASE_DIR + f'cube_Dfull_20260123_{seeing}_3600_100ms'
                          + f'_Kmag{scao_K}_0piston_meters_scao_only_{npupil}.fits')
    # For f_phase, the magnitude is floored to the minimum magnitude available in
    # the ALF sensor-noise files, so that bright stars share a single phase cube
    # with the faintest star for which ALF data exists.
    _min_alf_mag = get_min_alf_magnitude(band, mode, ncpa_freq,
                                         data_dir=dir_input+'/wavefront/alf/')
    fphase_mag   = max(magnitude, _min_alf_mag)
    f_phase  = (PHASE_DIR + f'cube_Dfull_20260123_{seeing}_3600_100ms'
                          + f'_Kmag{scao_K}_{band}mag{fphase_mag}_{mode}_all_{npupil}.fits')

    # PSF output directory: uses fphase_mag (not magnitude) because the PSF
    # is shared by all stars brighter than the ALF grid boundary.
    dir_output_psf = (f'output_files/2026_grid/'
                      f'{band}_{mode}_s={seeing}_mag={fphase_mag}'
                      f'_{duration}s_{dit*1e3:.0f}ms/')

    dir_output = (f'output_files/2026_grid/'
                  f'{band}_{mode}_s={seeing}_mag={magnitude}'
                  f'_{duration}s_{dit*1e3:.0f}ms/')

    dir_output_psf=os.path.join(dir_current, dir_output_psf)
    dir_output=os.path.join(dir_current, dir_output)

    # WFE noise scalars — computed based on NCPA parameters
    sigLF = None    # in nm rms 
    sigHF = None    # in nm rms 
    if do_f_phase:
        sigLF = get_wfe_tip_tilt(band) # in nm rms 
        sigHF = get_wfe_higher_order(band, magnitude, mode, ncpa_freq,
                                     nzern=ncpa['nmodes'],
                                     data_dir=dir_input+'/wavefront/alf/') # in nm rms 

def get_wfe_tip_tilt(band):
    """
    Return the tip-tilt WFE contribution from QACITS centering error [nm rms].

    The QACITS pointing error (stored in QACITS_ERROR_MAS, in mas peak-to-valley)
    is converted to nm rms using::

        wfe [nm rms] = error [mas]  *  (1e-3 / 206264.8)  *  D_TEL  *  ptv2rms

    where ``ptv2rms = 1/4`` converts peak-to-valley to rms, 
    and the angular-to-linear factor uses the ELT primary diameter.

    Parameters
    ----------
    band : str
        Observing band, one of 'L', 'M', 'N1', 'N2'.

    Returns
    -------
    float
        Tip-tilt WFE [nm rms] for the requested band.
    """    
    ptv2rms = 1 / 4
    mas2nm = 1e-3 / 206264.8 * D_TEL * ptv2rms * 1e9
    qacits_error_nm = {key: val * mas2nm for key, val in QACITS_ERROR_MAS.items()}

    # return sigLF
    return qacits_error_nm[band]


def get_min_alf_magnitude(band, mode, fr,
                 data_dir='/home/gorban/python_scripts/alf/results/'):
    """
    Return the minimum stellar magnitude available in the ALF sensor-noise
    FITS file for the given band, mask, and frame rate.

    This is used to floor the magnitude label in the phase-cube filename so
    that stars brighter than the ALF grid boundary share a single cached file.

    Parameters
    ----------
    band : str
    mask : str
        'CVC', 'RAVC'
    fr : float
        NCPA correction loop frame rate [Hz].
    data_dir : str, optional

    Returns
    -------
    float
        Minimum magnitude in the ALF table.
    """
    mask = MASK_TAGS[band][mode]
    sensor_noise_fname = data_dir + f'ALF_sensor_noise_{band}_{mask}_{fr}Hz.fits'
    mag_vec = fits.getdata(sensor_noise_fname)[0]
    return float(mag_vec[0])


def get_wfe_higher_order(band, magnitude, mode, fr,
                 nzern=20,
                 data_dir='/home/gorban/python_scripts/alf/results_data/',
                 min_clip=True):
    """
    Return the higher-order WFE contributed by ALF (NCPA) sensor noise [nm rms].

    Reads pre-computed ALF sensor-noise curves from FITS files and interpolates
    linearly to the requested stellar magnitude.  Each file contains a 2-row
    array: ``[mag_vector, sensor_noise_vector]``.

    Special cases
    -------------
    - RAVC / L-band: only valid for mag <= 4; assumed equal to the CLS-LM
      L-band residual for brighter stars.

    Parameters
    ----------
    band : str
        Observing band, one of 'L', 'M', 'N1', 'N2'.
    magnitude : float
        Stellar magnitude in the observing band.
    mode : str
        'CVC', 'RAVC'
    fr : float
        NCPA correction loop frame rate [Hz].
    nzern : int, optional
        Number of Zernike modes corrected by ALF (default 20).
        Does not include piston, tip, or tilt (as of 11/6/2026).
    data_dir : str, optional
        Directory containing the ALF sensor-noise FITS files.
    min_clip : bool, optional
        If True, clip the interpolated noise to the minimum tabulated value
        (prevents unphysical extrapolation at bright magnitudes).

    Returns
    -------
    float
        Higher-order WFE standard deviation per mode [nm rms].
    """
    mask = MASK_TAGS[band][mode]

    sensor_noise_fname = data_dir + f'ALF_sensor_noise_{band}_{mask}_{fr}Hz.fits'
    tmp = fits.getdata(sensor_noise_fname)

    mag_vec = tmp[0]
    sensor_noise = tmp[1]
    #   Find the two closest magnitudes
    idx = np.searchsorted(mag_vec, magnitude, side='right') - 1
    idx = np.clip(idx, 0, len(mag_vec) - 2) 
    # Perform linear interpolation
    x0, x1 = mag_vec[idx], mag_vec[idx + 1]
    y0, y1 = sensor_noise[idx], sensor_noise[idx + 1]

    slope = (y1 - y0) / (x1 - x0)
    interpolated_noise = y0 + slope * (magnitude - x0)

    if min_clip:
        # clippe noise to minimum from ALF simulation
        if interpolated_noise < sensor_noise[0]:
            interpolated_noise = sensor_noise[0]

    std_noise_per_mode = interpolated_noise / np.sqrt(nzern)
    # return sigHF
    return std_noise_per_mode


# --- Main Simulation Function ---
@ex.automain # immediately calls ex.run_commandline(), and then if __name__=='__main__'
# @ex.main
def run_simulation(mode, band, magnitude, duration, dit,
                   seeing, ncpa,
                   scao_K, wv_rms, ncpa_freq, sigLF, sigHF,
                   fphase_mag, dir_output, dir_output_psf,
                   f_phase, f_wv, f_cbw, f_scao, f_amp, f_oat, f_pupil,
                   do_f_phase, do_propagation, do_contrast_curves,
                   dry_run,
                   _config, _run):
    """
    Main Sacred experiment function for HEEPS contrast curve simulation.
    
    This function orchestrates a complete end-to-end HCI simulation including:
    1. Environment fingerprinting (HEEPS + VIP git info, package versions)
    2. Combined phase cube generation (SCAO + WV + NCPA with closed-loop correction)
    3. HEEPS wavefront propagation (pupil, on-axis and off-axis PSFs)
    4. Contrast curve computation (raw and post-processed with/without photon noise)
    5. Visualization and artifact logging to Sacred FileStorageObserver
    
    All parameters are provided automatically by Sacred from the experiment config
    (default_config, derived_config, and ingredient configs) and any command-line
    or programmatic overrides via config_updates.
    
    Parameters
    ----------
    mode : str
        HCI mode, one of 'CVC', 'RAVC', 'APP', 'SPP', 'CLC', 'ELT'.
    band : str
        Observing band, one of 'L', 'M', 'N1', 'N2'.
    magnitude : float
        Stellar magnitude in the observing band.
    duration : float
        Total ADI sequence duration [s].
    dit : float
        Detector integration time [s].
    seeing : str
        Atmospheric seeing percentile: 'Q1', 'Q2', 'Q3', or 'Q4'.
    ncpa : dict
        NCPA ingredient config with keys 'frequency' (Hz), 'lag' (frames),
        'gain_I' (float), 'nmodes' (int).
    scao_K : int
        SCAO K-magnitude (brightness-dependent, derived from magnitude and SCAO_K_THRES).
    wv_rms : float
        RMS water-vapour OPD [nm] (derived from band and seeing via WV_TABLE).
    ncpa_freq : float
        NCPA correction loop frequency [Hz] (derived from magnitude and NCPA_MAG_THRES).
    sigLF : float
        Tip-tilt WFE [nm rms] from QACITS centering error (derived in main config).
    sigHF : float
        Higher-order WFE per mode [nm rms] from ALF sensor noise (derived in main config).
    fphase_mag : float
        Magnitude label for phase cube filename (floored to min ALF magnitude).
    dir_output : str
        Output directory for this specific magnitude run.
    dir_output_psf : str
        Shared PSF output directory (magnitude-independent for caching).
    f_phase : str
        Combined phase cube filename (SCAO + WV + NCPA).
    f_wv : str
        Water-vapour phase cube filename.
    f_cbw : str
        CBW NCPA phase cube filename.
    f_scao : str
        SCAO residual phase cube filename.
    f_amp : str
        Talbot effect (amplitude screen) filename.
    f_oat : str
        Vortex off-axis transmission filename.
    f_pupil : str
        Entrance pupil mask filename.
    do_f_phase : bool
        If True, (re-)generate the combined phase cube if it doesn't exist.
    do_propagation : bool
        If True, run HEEPS wavefront propagation.
    do_contrast_curves : bool
        If True, compute contrast curves.
    dry_run : bool
        If True, skip all heavy computation (phase cube, propagation, contrast curves).
    _config : dict
        Sacred's complete configuration dictionary for this run.
    _run : sacred.run.Run
        Sacred Run object for adding artifacts and logging.
    
    Returns
    -------
    None
        All outputs are saved to disk and logged as Sacred artifacts:
        - fingerprint.txt: environment info
        - conf_{band}_{mode}.pkl: HEEPS configuration
        - cc_raw_{band}_{mode}.fits: raw contrast curve
        - cc_adi_bckg0_{band}_{mode}.fits: post-processed contrast (no photon noise)
        - cc_adi_bckg1_{band}_{mode}.fits: post-processed contrast (with photon noise)
        - contrast_curve.png: visualization (not logged to reduce redundancy)
    
    Notes
    -----
    - PSF files (pupil, onaxis_PSF, offaxis_PSF) are shared across magnitudes
      brighter than the ALF magnitude boundary and stored in dir_output_psf.
      Magnitude-specific runs create symlinks to these shared files.
    - If output files already exist, they are reused (cached) to avoid redundant
      computation.
    - The function is decorated with @ex.automain, so running the script directly
      executes one Sacred run with default/CLI config, then any __main__ block runs.
    
    See Also
    --------
    PhaseCubeGenerator : Generates combined phase cubes with NCPA correction.
    get_wfe_tip_tilt : Computes tip-tilt WFE from QACITS error.
    get_wfe_higher_order : Computes higher-order WFE from ALF sensor noise.
    """

    # print('*** TOTO DEBUG ***', flush=True)

    if dry_run:
        print('*** DRY-RUN MODE: skipping all heavy computation ***')
        do_f_phase = do_propagation = do_contrast_curves = False

    # Log the final configuration
    print(f"""
    Running simulation:
    - Mode: {mode}, Band: {band}, Magnitude: {magnitude}
    - Seeing: {seeing}
    - Water vapor: {wv_rms} nm rms
    - SCAO K-mag: {scao_K}
    - NCPA : lag={ncpa['lag']}, nmodes={ncpa['nmodes']}, {ncpa_freq} Hz
    - WFE noise: sigLF={sigLF:.1f} nm rms, sigHF={sigHF:.1f} nm rms, sigHF_all={np.sqrt(ncpa['nmodes'])*sigHF:.1f} nm rms
    -> f_phase : {f_phase}

    Input files:
    - f_pupil  : {f_pupil}
    - f_scao   : {f_scao}
    - f_wv     : {f_wv}
    - f_cbw    : {f_cbw}
    - f_amp    : {f_amp}
    - f_oat    : {f_oat}


    """)

    #*****************#
    # -- Fingerprint
    print(getFingerprint())

    tmp_fname = 'temp_fingerprint.txt'
    writeToFile('./', tmp_fname)
    # adding the file to the run artifact (stored by FileStorageObserver)
    _run.add_artifact(tmp_fname, name='fingerprint.txt')
    os.remove(tmp_fname)

    #*****************#
    print('\n *** GENERATE PHASE CUBE ***')
    print(f' f_phase = {f_phase}')
    if do_f_phase and not os.path.isfile(f_phase):
        generator = PhaseCubeGenerator(_config)
        generator.run(sigLF, sigHF)
    else:
        print('--- SKIPPING creation of f_phase ---')

    # ------------------------------------ #
    # Run HEEPS contrast curve calculation
    print('\n *** HEEPS simulation ***')
    print('---- CONFIGURATION ---')
    
    # Sanitize Sacred config for HEEPS (remove Sacred internal keys and non-serializable objects)
    heeps_config = {k: v for k, v in _config.items() 
                    if not k.startswith('_') and isinstance(v, (str, int, float, bool, type(None)))}
    
    # Map Sacred parameter names to HEEPS names
    if 'magnitude' in heeps_config:
        heeps_config['mag'] = heeps_config.pop('magnitude')
    
    conf = heeps.config.read_config(verbose=False, **heeps_config)
    conf = heeps.config.update_config(saveconf=True, verbose=True, **conf)
    _run.add_artifact(os.path.join(dir_output, f'conf_{band}_{mode}.pkl'), name='conf_heeps.pkl')

    _pupil_file   = os.path.join(dir_output_psf, f'pupil_{band}_{mode}.fits')
    _onaxis_file  = os.path.join(dir_output_psf, f'onaxis_PSF_{band}_{mode}.fits')
    _offaxis_file = os.path.join(dir_output_psf, f'offaxis_PSF_{band}_{mode}.fits')
    _psf_exists   = all(os.path.isfile(f) for f in [_pupil_file, _onaxis_file, _offaxis_file])

    print('---- PROPAGATION ---')
    if do_propagation:
        if not _psf_exists:
            print('-- Initialize pupil & wavefront --')
            conf['dir_output'] = dir_output_psf   # direct HEEPS output to shared PSF directory

            wf = heeps.pupil.pupil(savefits=True, verbose=True, **conf)

            # Create off-axis PSF
            heeps.wavefront.propagate(wf, onaxis=False, avg=True, savefits=True, verbose=True, **conf)

            # Create on-axis PSF --- computationally intensive !
            heeps.wavefront.propagate(wf, onaxis=True, savefits=True, verbose=True, **conf)

            conf['dir_output'] = dir_output       # reset to per-magnitude dir for contrast curves
        else:
            print('--- SKIPPING PROPAGATION (output files already exist) ---')

        print('Create symbolic links in dir_output pointing to the shared PSF files')
        os.makedirs(dir_output, exist_ok=True)
        for _src, _fname in [(_pupil_file,   f'pupil_{band}_{mode}.fits'),
                              (_onaxis_file,  f'onaxis_PSF_{band}_{mode}.fits'),
                              (_offaxis_file, f'offaxis_PSF_{band}_{mode}.fits')]:
            _link = os.path.join(dir_output, _fname)
            if not os.path.exists(_link):
                os.symlink(os.path.abspath(_src), _link)
                print(f'  symlink: {_link} -> {_src}')
    else:
        print('--- SKIPPING PROPAGATION (do_propagation = False) ---')

    if do_contrast_curves:
        print('---- CONTRAST CURVES ---')
        
        # Raw contrast
        cc_raw_file = os.path.join(dir_output, f'cc_raw_{band}_{mode}.fits')
        if os.path.isfile(cc_raw_file):
            print(' -- Raw contrast: loading existing file --')
            cc_raw_data = fits.getdata(cc_raw_file)
            sep = cc_raw_data[0]
            raw = cc_raw_data[1]
        else:
            print(' -- Raw contrast: computing --')
            sep, raw = heeps.contrast.cc_raw(savefits=True, verbose=True, **conf)
        _run.add_artifact(cc_raw_file, name='cc_raw.fits')

        # Post-processed 5-sigma contrast curve (without photon noise)
        conf['cpu_count'] = 1   # this is faster see email 12/3/2026
        print(f"  setting cpu_count={conf['cpu_count']} for cc_adi --")
        
        cc_adi_file = os.path.join(dir_output, f'cc_adi_bckg0_{band}_{mode}.fits')
        if os.path.isfile(cc_adi_file):
            print(' -- Post-processed contrast (w/o photon noise): loading existing file --')
            cc_adi_data = fits.getdata(cc_adi_file)
            sep1 = cc_adi_data[0]
            adi1 = cc_adi_data[1]
        else:
            print(' -- Post-processed contrast (w/o photon noise): computing --')
            sep1, adi1 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)
        _run.add_artifact(cc_adi_file, name='cc_adi.fits')

        # Post-processed contrast with star flux + background + photon noise
        conf['add_bckg'] = True
        cc_adi_bckg_file = os.path.join(dir_output, f'cc_adi_bckg1_mag{magnitude}_{band}_{mode}.fits')
        if os.path.isfile(cc_adi_bckg_file):
            print(' -- Post-processed contrast (with photon noise): loading existing file --')
            cc_adi_bckg_data = fits.getdata(cc_adi_bckg_file)
            sep2 = cc_adi_bckg_data[0]
            adi2 = cc_adi_bckg_data[1]
        else:
            print(' -- Post-processed contrast (with photon noise): computing --')
            sep2, adi2 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)
        _run.add_artifact(cc_adi_bckg_file, name='cc_adi_bckg.fits')

        # Create figure
        print('  -- Create figure -- ')
        savename = 'contrast_curve.png'
        xlabel = 'Angular separation $[arcsec]$'
        ylabel_adi = r'5-$\sigma$ sensitivity (contrast)'
        ylabel_raw = 'raw contrast'
        fig = plt.figure(figsize=(12, 6))
        fig.subplots(2, 1, sharex=True)
        fig.subplots_adjust(hspace=0.02)
        axes = fig.axes
        axes[0].set_ylim(1e-6, 1e-2)
        axes[1].set_ylim(1e-8, 1e-3)
        axes[1].set_xlabel(xlabel)
        for j, (ax, ylabel) in enumerate(zip(axes, [ylabel_raw, ylabel_adi])):
            ax.set_ylabel(ylabel)
            ax.grid(True), ax.grid(which='minor', linestyle=':')
            ax.loglog()
            ax.xaxis.set_major_formatter(plt.ScalarFormatter())
            # ax.set_xticks([0.06, 0.1, 0.2, 0.5, 1, 1.2])  # N-band
            # ax.set_xlim(0.06, 1.2)
            if band=='N2' or band=='N1':
                ax.set_xticks([0.06, 0.1, 0.2, 0.5, 1, 1.2]) # N-band
                ax.set_xlim(0.06, 1.2)
            else:
                ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
                ax.set_xlim(0.02, 0.8)

        axes[0].plot(sep, raw, 'C0', label=mode, marker='d', markevery=0.12, markersize=4)
        axes[0].legend()
        axes[0].set_title(f'{mode} {band}-mag = {magnitude}')
        axes[1].plot(sep1, adi1, 'C0', label=mode, marker='d', markevery=0.12, markersize=4)
        axes[1].plot(sep2, adi2, ':C0', label='background', marker='d', markevery=0.12, markersize=4)
        axes[1].legend(ncol=2, loc='upper right')
        fig.savefig(f"{conf['dir_output']}/{savename}", dpi=300, transparent=True)
        # -- removing png from db -> redundant with the cc_adi files that are already copy to db
        # _run.add_artifact(f"{conf['dir_output']}/{savename}", name='contrast_curve.png')


# if __name__=="__main__":
#     from multiprocessing import get_context
#     # import subprocess
#     import sys
#     from itertools import product
#     import numpy as np

#     # Define magnitude ranges per band
#     dmag=0.5
#     band_magnitude_ranges = {
#         # "L": np.arange(-1.5, 9+dmag, dmag),  # Example: L-band magnitudes
#         # "M": np.arange(-1.5, 7+dmag, dmag),       # Example: M-band magnitudes
#         # "N1": np.arange(-1.5, 3+dmag, dmag),        # Example: N-band magnitudes
#         # "N2": np.arange(-1.5, 3+dmag, dmag)        # Example: N-band magnitudes
#         "L": [6],  # Example: L-band magnitudes
#     }

#     seeings = ["Q1", "Q2", "Q3"] # ESO definition of percentile (median is between Q2 and Q3)
#     seeings = ["Q2"] # ESO definition of percentile (median is between Q2 and Q3)


#     # Generate all combinations dynamically
#     grid_configs = []
#     for band, magnitudes in band_magnitude_ranges.items():
#         for magnitude, seeing in product(magnitudes, seeings):
#             grid_configs.append({
#                 "band": band,
#                 "magnitude": magnitude,
#                 "seeing": seeing,
#                 "dry_run":False
#             })

#     # # Run each configuration sequentially with stdout capture
#     for config in grid_configs:
#         ex.run(config_updates=config, options={'--capture': 'sys'})

#     # # Function to run a single configuration
#     # def run_single_config(config):
#     #     result = ex.run(config_updates=config)

#     # # Run in parallel with 4 jobs --> DOES NOT WORK, because create_phase_map_all.py uses multiprocesses !!
#     # n_job_max = 4
#     # with get_context('fork').Pool(n_job_max) as pool:
#     #     pool.map(run_single_config, grid_configs)



