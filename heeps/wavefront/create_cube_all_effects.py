from heeps.pupil.create_petal import create_petal
from heeps.util.freq_decomp import conv_kernel, spatial, temporal
from heeps.config.read_config import read_config
from heeps.config.update_config import update_config
from heeps.util.img_processing import resize_cube
from astropy.io import fits
import numpy as np
import os

# inputs
os.chdir(os.path.normpath(os.path.expandvars('$HOME/heeps_metis/input_files')))
nimg = 511
npupil = 285 # L band
cpp = 10
pupil_img_size = 39.9988
diam_nominal = 38.542

master_seed = {'LTF':123456, 'HTF':345678}
cutoff = 0.01 # in Hz
tag = 'Cbasic_20201130'
t_max = 3600 # in s
dt = 300 # in ms
npetals = 6
f_scao_screens = 'wavefront/cbasic/cube_%s_%ss_%sms_0piston_meters_%s_%s_%s.fits'
f_mask = 'wavefront/cbasic/mask_%s_%s_%s.fits'
f_ncpa_frame = 'wavefront/ncpa/DIFF_LM_20201119_fullM1.fits'
f_ncpa_screens = 'wavefront/ncpa/ncpa_fullM1_%s_%scpp_%s_%s.fits'
f_petal_screens = 'wavefront/ncpa/cube_petal_piston_%s_seed=%s_%s_%s.fits'

# calculate cube size
nframes = int(t_max/dt*1000)
# load scao atmosphere residuals and mask
scao = fits.getdata(f_scao_screens%(tag, t_max, dt, 'scao_only', 'L', npupil))
mask = fits.getdata(f_mask%(tag, 'L', npupil))

for rep in [1]:#np.arange(1,11):
    print('rep=%s'%rep, end=', ')

    # create ncpa maps (spatial frequencies)
    try:
        ncpa_allSF = fits.getdata(f_ncpa_screens%('allSF', cpp, 'L', npupil))
        ncpa_LSF = fits.getdata(f_ncpa_screens%('LSF', cpp, 'L', npupil))
        ncpa_HSF = fits.getdata(f_ncpa_screens%('HSF', cpp, 'L', npupil))
        print('NCPA spatial files loaded.')
    except FileNotFoundError:
        print('Creating NCPA spatial files...')
        ncpa_frame = fits.getdata(f_ncpa_frame)
        nkernel = nimg*diam_nominal/pupil_img_size
        kernel = conv_kernel(nkernel, cpp)
        ncpa_allSF, ncpa_LSF, ncpa_HSF = spatial(ncpa_frame, kernel, npupil=npupil, norm=True, verbose=True)
        fits.writeto(f_ncpa_screens%('allSF', cpp, npupil), np.float32(ncpa_allSF))
        fits.writeto(f_ncpa_screens%('LSF', cpp, npupil), np.float32(ncpa_LSF))
        fits.writeto(f_ncpa_screens%('HSF', cpp, npupil), np.float32(ncpa_HSF))
    finally:
        LTF1 = fits.getdata('wavefront/ncpa/time_series_LTF_0-0.01Hz_12000x1rms_seed=832404.fits')
        LTF2 = fits.getdata('wavefront/ncpa/time_series_LTF_0-0.01Hz_12000x1rms_seed=523364.fits')
        HTF1 = fits.getdata('wavefront/ncpa/time_series_HTF_0.01-1Hz_12000x1rms_seed=832404.fits')
        ncpa_STA = np.array([ncpa_HSF]*nframes)
        ncpa_QLSF = np.array([x*ncpa_LSF for x in LTF1])
        ncpa_QHSF = np.array([x*ncpa_HSF for x in LTF2])
        ncpa_DYN = np.array([x*ncpa_allSF for x in HTF1])

    # create petal piston cubes
    try:
        f_piston_Q = 'wavefront/ncpa/cube_petal_piston_%s_seed=%s_%s_%s.fits'\
            %('LTF', master_seed['LTF'], 'L', npupil)
        f_piston_DYN = 'wavefront/ncpa/cube_petal_piston_%s_seed=%s_%s_%s.fits'\
            %('HTF', master_seed['HTF'], 'L', npupil)
        piston_Q = fits.getdata(f_piston_Q)
        piston_DYN = fits.getdata(f_piston_DYN)
        print('Petal piston files loaded.')
    except FileNotFoundError:
        print('Creating petal piston files...')
        petals = np.float32([create_petal(x, npetals, npupil) for x in range(npetals)])
        for TF, fc1, fc2 in zip(['LTF', 'HTF'], [0, cutoff], [cutoff, 1]):
            # create time series with random seeds
            np.random.seed(master_seed[TF])
            seeds = np.random.randint(1e5, 1e6, 6)
            tseries = np.float32([temporal(t_max, dt/1000, fc1, fc2, seed=seed) for seed in seeds])
            # print rms 
            pp_rms = np.mean(np.std(tseries, 0))
            print('mean %s = %3.2f nm rms'%(TF, pp_rms))
            # multiply time series by their respective petal, and sum them up
            pp = np.sum([np.array(petals[x], ndmin=3).T * tseries[x] \
                for x in range(npetals)], 0).T
            # mask, normalize and save petal piston
            pp *= mask
            pp /= pp_rms
            fits.writeto(f_petal_screens%(TF, master_seed[TF], 'L', npupil), np.float32(pp), overwrite=True)
        piston_Q = fits.getdata(f_piston_Q)
        piston_DYN = fits.getdata(f_piston_DYN)

    # All effects
    ncpa_QLSF = (ncpa_QLSF + piston_Q)/np.sqrt(2)
    ncpa_QHSF = ncpa_QHSF 
    ncpa_DYN = (ncpa_DYN + piston_DYN)/np.sqrt(2)
    # norm (almost 1)
    ncpa_QLSF /= np.mean([np.nanstd(x) for x in ncpa_QLSF])
    ncpa_QHSF /= np.mean([np.nanstd(x) for x in ncpa_QHSF])
    ncpa_DYN /= np.mean([np.nanstd(x) for x in ncpa_DYN])
    # total in meters (values at L band)
    ncpa_piston_ALL = ncpa_STA*36e-9 + ncpa_QLSF*20e-9 + ncpa_QHSF*20e-9 + ncpa_DYN*40e-9
   # ncpa_piston_ALL = ncpa_QLSF*20e-9 + ncpa_QHSF*20e-9 + ncpa_DYN*40e-9
   # ncpa_piston_ALL = ncpa_QLSF*20e-9 + ncpa_QHSF*20e-9

    # TOTAL PHASE SCREENS
    lamL = update_config(**dict(read_config(), band='L'))['lam']
    for band in ['L', 'M', 'N1', 'N2']:
        conf = update_config(**dict(read_config(), band=band))
        scaling = conf['lam']/lamL
        #print(scaling, conf['npupil'])
        f_out = f_scao_screens%(tag, t_max, dt, 'all_ncpa', band, conf['npupil'])
        fits.writeto(f_out, resize_cube(scao + ncpa_piston_ALL*scaling, conf['npupil']), overwrite=True)
        print('%s created'%f_out)