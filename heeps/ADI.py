#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import vip_hci
import numpy as np
from astropy.io import fits
import os.path


""" user inputs """
folder = '/Volumes/Data/METIS/compass600s10ms/'
#folder = '/mnt/disk4tb/METIS/compass600s10ms/'
mode = 'RAVC'                   # modes: ELT, VC, RAVC, APP, CL_4, CL_5...
scao_name = 'compass'           # SCAO simulator name
cube_duration = 600             # SCAO cube duration in seconds
cube_samp = 100                 # SCAO cube sampling in ms
adi_cube_duration = 3600        # ADI cube duration in seconds
adi_cube_samp = 100             # ADI cube sampling in ms
adi_cube_avg = 0                # ADI cube averaging in ms
lat = -24.63                    # telescope latitude in deg (Paranal -24.63)
dec = -2.47                     # star declination in deg (e.g. 51 Eri -2.47)
Lmag = 5                        # star magnitude at L
mag5_ADU = 3.4539e9             # L'=5 star flux in counts(ADU) per sec (@NACO L' filter)
bckg_ADU = 229.78e3             # background flux in counts(ADU) per sec (@NACO L' filter)
rim = 19                        # psf image radius (in pixels)
psc_simu = 0.005                # simulations (SCAO) platescale in arcsec/pix
psc_inst = 0.00525              # instruments (METIS) platescale in arcsec/pix
calc_trans = False              # true if transmission must be calculated
plot_cc = True                  # true if plot contrast curve
algo = vip_hci.medsub.median_sub# VIP post-processing algorithm

# output filename (can carry some values)
filename = '%s%ss_samp%sms_ADI%ss_samp%sms_avg%sms_dec%sdeg_magL%s_%s' \
        %(scao_name, cube_duration, cube_samp, adi_cube_duration, adi_cube_samp, \
        adi_cube_avg, dec, Lmag, mode)
#filename = 'testname'

""" transmission : ratio of intensities (squared amplitudes) in Lyot-Stop plane """
if calc_trans is True:
    I_ELT = fits.getdata(os.path.join(folder, 'ELT_LS.fits'))**2
    I_OFFAXIS = fits.getdata(os.path.join(folder, 'OFFAXIS_%s_LS.fits'%mode))**2
    trans = np.sum(I_OFFAXIS)/np.sum(I_ELT)
else:
    trans_all = {'ELT': 1.,
                  'VC': 0.9012406091763115,
                'RAVC': 0.3392759549914341,
                'CL_4': 0.502257180317745} # 4-lam/D diam classical Lyot mask
    trans = trans_all[mode]

""" normalize off-axis and on-axis (cube) psfs """
# total flux of the non-coronagraphic PSF
psf_ELT = fits.getdata(os.path.join(folder, 'ELT_PSF.fits'))
ELT_flux = np.sum(psf_ELT)
# normalized off-axis PSF
psf_OFF = fits.getdata(os.path.join(folder, 'OFFAXIS_%s_PSF.fits'%mode))
psf_OFF /= ELT_flux
# normalized coronagraphic (on-axis) PSFs cube
psf_ON = fits.getdata(os.path.join(folder, 'ONAXIS_%s_PSF.fits'%mode))
psf_ON /= ELT_flux

""" VIP: resample psfs to match the instrument pixelscale """
psf_ON = vip_hci.preproc.cube_px_resampling(psf_ON, psc_simu/psc_inst)
psf_OFF = vip_hci.preproc.frame_px_resampling(psf_OFF, psc_simu/psc_inst)
# clip negative values due to resampling
psf_ON = psf_ON.clip(min=0)
psf_OFF = psf_OFF.clip(min=0)

""" add noises: star, background, photon """
# cube shape
(ncube, xon, yon) = psf_ON.shape
# detector integration time
DIT = adi_cube_duration/ncube
# star flux
star_flux = DIT * mag5_ADU * 10**(-0.4*(Lmag - 5))
# background flux
bckg_flux = DIT * bckg_ADU
# add the stellar flux and the background
psf_OFF *= star_flux
psf_ON *= star_flux
psf_ON += bckg_flux*trans
# generate a cube of random noise ~ N(0,1) * sqrt(psf)
noise = np.random.randn(ncube, xon, yon) * np.sqrt(psf_ON)
psf_ON += noise

""" VIP: aperture photometry of psf_OFF used to scale the contrast """
# get the center pixel
(xoff, yoff) = psf_OFF.shape
(cx, cy) = (int(xoff/2), int(yoff/2))
# fit a 2D Gaussian --> output: fwhm, x-y centroid
fit = vip_hci.var.fit_2dgaussian(psf_OFF[cx-rim:cx+rim+1, cy-rim:cy+rim+1], True, \
        (rim,rim), debug=False, full_output=True)
# derive the FWHM
fwhm = np.mean([fit['fwhm_x'],fit['fwhm_y']])
# recenter and crop
psf_OFF = vip_hci.preproc.frame_shift(psf_OFF, rim-fit['centroid_x'], rim-fit['centroid_y'])
psf_OFF_crop = psf_OFF[cx-rim:cx+rim+1, cy-rim:cy+rim+1]
# FWHM aperture photometry of psf_OFF_crop
starphot = vip_hci.metrics.aperture_flux(psf_OFF_crop,[rim],[rim],fwhm)[0]

""" parallactic angles for ADI """
# duration -> hour angle conversion
ha = adi_cube_duration/3600/24*360
# angles in rad
hr = np.deg2rad(np.linspace(-ha/2, ha/2, ncube))
dr = np.deg2rad(dec)
lr = np.deg2rad(lat)
# parallactic angle in deg
pa = -np.rad2deg(np.arctan2(-np.sin(hr), np.cos(dr)*np.tan(lr)-np.sin(dr)*np.cos(hr)))

""" VIP: post-processing (ADI, ADI-PCA,...) """
# psf after post-processing
out, derot, psf_pp = algo(psf_ON, pa, full_output=True)
# contrast curve after post-processing
cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop, fwhm, psc_inst, \
        starphot, algo=algo, nbranch=1, sigma=5, debug=False, plot=False)

""" saving to fits files """
fits.writeto(os.path.join(folder, 'psf_' + filename + '.fits'), psf_pp, overwrite=True)
hdu = fits.PrimaryHDU(cc_pp)
hdu.writeto(os.path.join(folder, 'cc_' + filename + '.fits'), overwrite=True)

""" figure """
if plot_cc is True:
    x = cc_pp.get('distance_arcsec')
    y1 = cc_pp.get('sensitivity_gaussian')
    y2 = cc_pp.get('sensitivity_student')
    plt.figure(1)
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.yscale('log')
    plt.grid('on')
    plt.xlabel("Angular separation [arcsec]")
    plt.ylabel(r"5-$\sigma$ sensitivity")
    plt.legend()
    plt.xlim([0, 0.6])
    plt.ylim([1e-7, 1e-0])
    plt.show(block=False)
    plt.savefig(os.path.join(folder, 'cc_' + filename + '.png'), dpi=300, transparent=True)
