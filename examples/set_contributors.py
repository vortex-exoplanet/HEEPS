import numpy as np
import astropy.units as u
from astropy.io import fits
from heeps.config import conf
import os

#----------------------------------#
#  ALL CONTRIBUTORS
conf['input_dir'] = '/mnt/disk4tb/METIS/heeps_analysis/input_files'

# load ncpa maps (spatial frequencies)
ncpa_allSF = fits.getdata(os.path.join(conf['input_dir'],\
                                       'ncpa_allSF_253.fits'))
ncpa_LSF = fits.getdata(os.path.join(conf['input_dir'],\
                                     'ncpa_LSF_10cpp_253.fits'))
ncpa_HSF = fits.getdata(os.path.join(conf['input_dir'],\
                                     'ncpa_HSF_10cpp_253.fits'))
mean_allSF = np.nanmean(ncpa_allSF)
mean_LSF = np.nanmean(ncpa_LSF)
mean_HSF = np.nanmean(ncpa_HSF)

# load time series (temporal frequencies)
LTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=832404.fits'))
LTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=523364.fits'))
LTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=409566.fits'))
LTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=224788.fits'))
HTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=832404.fits'))
HTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=523364.fits'))
HTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=409566.fits'))
HTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=224788.fits'))

# * NCPA cubes *
#   - STA : Static
#   - QLSF: quasi-static low spatial frequencies
#   - QHFS: quasi-static high spatial frequencies
#   - DYN : dynamic, ie high temporal frequencies
ncpa_STA = np.array([ncpa_HSF]*12000)
ncpa_QLSF = np.array([x*ncpa_LSF for x in LTF1])
ncpa_QHSF = np.array([x*ncpa_HSF for x in LTF2])
ncpa_DYN = np.array([x*ncpa_allSF for x in HTF1])
# Total NCPA with normalization
ncpa_ALL = ncpa_STA*35.9 + ncpa_QLSF*20 + ncpa_QHSF*20 + ncpa_DYN*40

# * Petal piston *
piston_QLSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_LSF_1rms_seed=123456.fits'))
piston_QHSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_HSF_1rms_seed=234567.fits'))
piston_DYN = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_HTF_allSF_1rms_seed=345678.fits'))

# * ncpa + piston *
ncpa_piston_QLSF = (ncpa_QLSF + piston_QLSF)/np.sqrt(2)  # %sqrt(2) to get 1nm rms
ncpa_piston_QHSF = (ncpa_QHSF + piston_QHSF)/np.sqrt(2)
ncpa_piston_DYN  = (ncpa_DYN + piston_DYN)/np.sqrt(2)
# norm (almost 1)
ncpa_piston_QLSF /= np.mean([np.nanstd(x) for x in ncpa_piston_QLSF])
ncpa_piston_QHSF /= np.mean([np.nanstd(x) for x in ncpa_piston_QHSF])
ncpa_piston_DYN  /= np.mean([np.nanstd(x) for x in ncpa_piston_DYN])
# total
ncpa_piston_ALL = ncpa_STA*35.9 + ncpa_piston_QLSF*20 +\
                ncpa_piston_QHSF*20 + ncpa_piston_DYN*40


# * Apodizer drift (1% ptv) *
#  pupil_drift = [[x,0,0,0,0,0] for x in np.linspace(-0.005,0.005,nframes)]
# * Apodizer drift (3% ptv) -> new baseline *
nframes=12000
pupil_drift = [[x,0,0,0,0,0] for x in np.linspace(-0.015,0.015,nframes)]

# * Pointing errors (zernikes [2,3]) *
#    - point_drift_x = np.linspace(-0.2, 0.2, nframes)
#    - point_jit_x = (np.random.normal(0, 2,\
#                      nframes)*u.mas).to('rad').value/(3.8e-6/37)
#      (GOX) :  conversion correct only for L-band !

mas2lamD = (1*u.mas).to('rad').value/(3.8e-6/37)    # mas -> lam/D
point_QSTA = np.array([LTF3, LTF4]).T * mas2lamD / np.sqrt(2) # in xy
point_DYN = np.array([HTF3, HTF4]).T * mas2lamD / np.sqrt(2)  # in xy
point_ALL = point_QSTA * 0.4 + point_DYN * 2       # factors in mas
