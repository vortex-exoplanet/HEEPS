#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import vip_hci
import numpy as np
from astropy.io import fits
import os.path


""" user inputs """
folder = '/Volumes/Data/METIS/compass600s10ms/'
#folder = '/mnt/disk4tb/METIS/compass600s10ms/'
mode = 'CL_4'                   # modes: ELT, VC, RAVC, APP, CL_4, CL_5...
scao_name = 'compass'           # SCAO simulator name
cube_duration = 600             # SCAO cube duration in seconds
cube_samp = 10                  # SCAO cube sampling in ms
adi_cube_duration = 3600        # ADI cube duration in seconds
adi_cube_samp = 10              # ADI cube sampling in ms
adi_cube_avg = 0                # ADI cube averaging in ms
dparall = 40                    # parallactic angle variation in deg
Lmag = 5                        # star magnitude here
mag5_ADU = 3.4539e9             # L'=5 star flux in counts(ADU) per sec (@NACO L' filter)
bckg_ADU = 229.78e3             # background flux in counts(ADU) per sec (@NACO L' filter)
rim = 19                        # psf image radius (in pixels)
psc_simu = 0.005                # simulations (SCAO) platescale in arcsec/pix
psc_inst = 0.00525              # instruments (METIS) platescale in arcsec/pix
calc_trans = False              # true if transmission must be calculated
plot_cc = True                  # true if plot contrast curve

# output filename (can carry some values)
filename = '%s%ss_samp%sms_ADI%ss_samp%sms_avg%sms_par%sdeg_magL%s_%s' \
        %(scao_name, cube_duration, cube_samp, adi_cube_duration, adi_cube_samp, \
        adi_cube_avg, dparall, Lmag, mode)
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

""" add noises: star, background, photon """
# cube length
ncube, x, y = psf_ON.shape
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
(x,y,z) = psf_ON.shape
noise = np.random.randn(x,y,z) * np.sqrt(psf_ON)
psf_ON += noise

""" VIP: resample psf_ON and psf_OFF to match the instrument pixelscale """
psf_ON = vip_hci.preproc.cube_px_resampling(psf_ON, psc_simu/psc_inst)
psf_OFF = vip_hci.preproc.frame_px_resampling(psf_OFF, psc_simu/psc_inst)
# clip negative values due to resampling
psf_ON = psf_ON.clip(min=0)
psf_OFF = psf_OFF.clip(min=0)

""" VIP: aperture photometry of psf_OFF used to scale the contrast """
# get the center pixel
(nx, ny) = psf_OFF.shape
(cx, cy) = (int(nx/2), int(ny/2))
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

""" VIP: ADI, ADI-PCA """
# parallactic angles during observation
angs = np.linspace(-dparall/2., dparall/2., ncube)
cc_adi = vip_hci.metrics.contrast_curve(psf_ON[:10], angs[:10], psf_OFF_crop, \
        fwhm, psc_inst, starphot, vip_hci.medsub.median_sub, nbranch=1, sigma=5, \
        debug=False, plot=False)
#out, derot, psf_ON_adi = vip_hci.medsub.median_sub(psf_ON, angs, full_output=True) # ADI
#psf_ON_pca = vip_hci.pca.pca(psf_ON, angs, ncomp=10, svd_mode='randsvd', \
#        full_output=False, check_mem=False) # ADI-PCA

""" saving to fits files """
hdu = fits.PrimaryHDU(cc_adi)
hdu.writeto(os.path.join(folder, 'cc_' + filename + '.fits'), overwrite=True)
#fits.writeto(os.path.join(folder, 'psf_' + filename + '.fits'), psf_ON_adi, overwrite=True)

""" figure """
if plot_cc is True:
    x = cc_adi.get('distance_arcsec')
    y1 = cc_adi.get('sensitivity_gaussian')
    y2 = cc_adi.get('sensitivity_student')
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
