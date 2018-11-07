#!/usr/bin/env python3

import vip_hci
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


out_dir = str('./output_files/')

#### inputs psfs and cube
psf_noMask = vip_hci.fits.open_fits(out_dir + 'test_PSF_ELT.fits')    # non-coronagraphic PSF for normalization (just a ELT PSF)
psf = vip_hci.fits.open_fits(out_dir +'test_PSF_OFFAXIS.fits')       # off-axis PSF (non-normalized) (No mask--with apodizer and LS)
#cube = vip_hci.fits.open_fits(out_dir + 'Heidelberg_100_PSF_cube_VC.fits')          # coronagraphic PSFs cube (non-normalized)
#cube = vip_hci.fits.open_fits(out_dir + 'Heidelberg_1000_PSF_cube_VC.fits')         # coronagraphic PSFs cube (non-normalized)
#cube = vip_hci.fits.open_fits(out_dir + 'Linz_100_PSF_cube_VC.fits')          # coronagraphic PSFs cube (non-normalized)
cube = vip_hci.fits.open_fits(out_dir + 'test_PSF_cube_VC.fits')[0:10]          # coronagraphic PSFs cube (non-normalized)


length,x,y = cube.shape

psf= psf/np.sum(psf_noMask)                                                 # off-axis PSF normalized wrt the total flux of the non-coronagraphic PSF
cube = cube/np.sum(psf_noMask)                                              # cube normalized wrt the total flux of the non-coronagraphic PSF

angs = np.linspace(-20, 20, length)
#angs = (np.linspace(160, 200, length) + 180) % 360 - 180

exp_t = 60*60/length                                                            # exposure time in sec for 1 hr sequence
Lmag = 5                                                                       # choose the star magnitude here
star_ref = exp_t * 3.4539e9                                                    # flux in ADU corresponding to a 5th mag star in the NACO-Lprime filter (METIS with no coronagraph)
star_flux = star_ref * 10**(-0.4*(Lmag-5))                                     # actual star flux for the chosen magnitude
bckg = exp_t * 229.78e3                                                        # flux in ADU/pix corresponding to the L-band background in the NACO-Lprime filter

#### pixelscale (arcsec/pix)
psc_simu = 0.005 # plate scale (arcsec/pix) in the simulations
psc_metis = 0.00525 # plate scale (arcsec/pix) of the METIS L-band detector



# to calculate the transmission
trans = 0.34  # if VC then use '0.81' and '0.34' in case of RAVC 

# Resample the disk PSF and cube with the pixelscale of the Metis instrument
cube_metis = vip_hci.preproc.cube_px_resampling(cube, psc_simu/psc_metis)
cube_metis = cube_metis.clip(min=0) # to avoid negative values due to the resampling
#cube_metis[cube_metis < 0] = 0 # to avoid negative values due to the resampling
psf_metis = vip_hci.preproc.frame_px_resampling(psf, psc_simu/psc_metis)


xpsf = int((psf_metis.shape[0])/2)
ypsf = int((psf_metis.shape[0])/2)

w = 19
fit = vip_hci.var.fit_2dgaussian(psf_metis[xpsf-w:xpsf+w+1, ypsf-w:ypsf+w+1], True, (w,w), debug=False, full_output=True)

# cropping the off-axis PSF
fwhmx = fit['fwhm_x']
fwhmy = fit['fwhm_y']
xcen = fit['centroid_x']
ycen = fit['centroid_y']

psf_cen = vip_hci.preproc.frame_shift(psf_metis, w-xcen, w-ycen)
psf_crop = psf_cen[xpsf-w:xpsf+w+1, ypsf-w:ypsf+w+1]

# Derive the FWHM
fwhm = np.mean([fwhmx,fwhmy])

# Adding the stellar flux and the background
psf_phot = psf_crop * star_flux
cube_phot = cube_metis * star_flux
cube_phot = cube_phot + bckg*trans # need to multiply the backgound fo the transmission

ran = np.random.randn(cube_phot.shape[0],cube_phot.shape[1],cube_phot.shape[2])
cube_phot = cube_phot + ran * np.sqrt(cube_phot)

starphot = vip_hci.metrics.aperture_flux(psf_phot,[w],[w],fwhm)

# ADI, ADI-PCA
cube_out_adi_mean100ms, cube_der_adi_mean100ms, image_adi_mean100ms = vip_hci.medsub.median_sub(cube_phot, angs, full_output=True) # ADI

ravc2_mean100ms_cc_adi = vip_hci.metrics.contrast_curve(cube_phot, angs, psf_crop, fwhm, psc_metis, starphot[0], vip_hci.medsub.median_sub, nbranch=1, debug=False, plot=False) # ADI contrast curves

x = ravc2_mean100ms_cc_adi.get('distance_arcsec')
y1 = ravc2_mean100ms_cc_adi.get('sensitivity_gaussian')
y2 = ravc2_mean100ms_cc_adi.get('sensitivity_student')

#np.savetxt('Heidelberg_1000_y1.txt', y1, delimiter=',')
#np.savetxt('Heidelberg_1000_x.txt', x, delimiter=',')

#np.savetxt('Linz_1000_y1.txt', y1, delimiter=',')
#np.savetxt('Linz_1000_x.txt', x, delimiter=',')

#
plt.figure(1, figsize=(12,8))
ax = plt.axes()
ax.semilogy(x, y1)
ax.semilogy(x,y2)
ax.grid()
ax.legend()
ax.set_ylabel("Contrast")
ax.set_xlabel("Arcsec")
ax.set_xlim([0, 0.6])
plt.show(block=False)
plt.savefig(out_dir + 'ADI.png')