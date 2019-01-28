import matplotlib.pyplot as plt
import vip_hci
import numpy as np
from astropy.io import fits
import os.path


def adi(path_offaxis='offaxis', path_onaxis='onaxis', path_output='output_files', 
        mode='VC', scao_name='compass', cube_duration=600, cube_samp=100, 
        adi_cube_duration=3600, adi_cube_samp=100, adi_cube_avg=0, lat=-24.63, 
        dec=-2.47, band='Lp', mag=5, psc_simu=5.21, psc_inst=5.21, rim=19, 
        add_bckg=True, calc_trans=False, plot_cc=False):
    """ 
    This function calculates and draws the contrast curve (5-sigma sensitivity) 
    for a specific set of off-axis PSF and on-axis PSF (or cube of PSFs),
    using the VIP package to perform ADI post-processing.
    
    Args:
        path_offaxis (string):
            Path to off-axis PSFs
        path_onaxis (string):
            Path to on-axis (cubes of) PSFs
        path_output (string):
            Path to output files
        mode (string):
            HCI mode: ELT, VC, RAVC, APP, CL4, CL5,...
        scao_name (string):
            SCAO simulator name
        cube_duration (integer):
            SCAO cube duration in seconds
        cube_samp (integer):
            SCAO cube sampling in ms
        adi_cube_duration (integer):
            ADI cube duration in seconds
        adi_cube_samp (integer):
            ADI cube sampling in ms
        adi_cube_avg (integer):
            ADI cube averaging in ms
        lat (float):
            Telescope latitude in deg (Paranal -24.63)
        dec (float):
            Star declination in deg (e.g. 51 Eri -2.47)
        band ():
            Spectral band ('Lp', 'Mp', 'N1', 'N2')
        mag (float):
            Star magnitude at band
        psc_simu (float):
            Simulation (fits files) pixel scale in mas/pix
        psc_inst (float):
            Instrument pixel scale in mas/pix (e.g. METIS LM=5.21, NQ=10.78)
        rim (integer):
            Psf image radius in pixels
        add_bckg (boolean)
            True means background flux and photon noise are added 
        calc_trans (boolean):
            True means transmission is calculated from fits files
        plot_cc (boolean):
            True means contrast curve is plotted
    """
    
    """ output filename """
    filename = '%s%ss_samp%sms_ADI%ss_samp%sms_avg%sms_dec%sdeg_%s_mag%s_bckg%s_%s' \
            %(scao_name, cube_duration, cube_samp, adi_cube_duration, \
            adi_cube_samp, adi_cube_avg, dec, band, mag, int(add_bckg), mode)
    
    """ transmission : ratio of intensities (squared amplitudes) in Lyot-Stop plane """
    if calc_trans is True:
        I_ELT = fits.getdata(os.path.join(path_offaxis, 'LS_%s_%s.fits' \
                %('ELT', band)))**2
        I_OFFAXIS = fits.getdata(os.path.join(path_offaxis, 'LS_%s_%s.fits' \
                %(mode, band)))**2
        trans = np.sum(I_OFFAXIS)/np.sum(I_ELT)
    else:
        trans_all = {'ELT': 1.,
                      'VC': 0.9012406091763115,
                    'RAVC': 0.3392759549914341,
                     'CL4': 0.502257180317745} # 4-lam/D diam classical Lyot mask
        trans = trans_all[mode]
    
    """ get normalized off-axis PSF (single) """
    # total flux of the non-coronagraphic PSF
    psf_ELT = fits.getdata(os.path.join(path_offaxis, 'PSF_%s_%s.fits' \
            %('ELT', band)))
    ELT_flux = np.sum(psf_ELT)
    # normalized off-axis PSF
    psf_OFF = fits.getdata(os.path.join(path_offaxis, 'PSF_%s_%s.fits' \
            %(mode, band)))
    psf_OFF /= ELT_flux
    
    """ get normalized on-axis PSFs (cube, resampled, and averaged) """
    # load cube, and format to 3D
    psf_ON = fits.getdata(os.path.join(path_onaxis, 'PSF_%s_%s.fits' \
            %(mode, band)))
    if psf_ON.ndim != 3:
        psf_ON = np.array(psf_ON, ndmin=3)
    # save PSF initial shape
    (ncube, xon, yon) = psf_ON.shape
    # resample based on ADI sampling vs simulation sampling
    nsamp = int(adi_cube_samp/cube_samp + .5)
    psf_ON = psf_ON[::nsamp,:]
    # average frames
    navg = max(1, int(adi_cube_avg/adi_cube_samp + .5))
    if navg > 1:
        end = ncube - ncube % navg
        psf_ON = np.mean(psf_ON[:end].reshape(-1, navg, xon, yon), axis=1)
    # normalized coronagraphic (on-axis) PSFs
    psf_ON /= ELT_flux
    
    """ VIP: resample psfs to match the instrument pixelscale """
    psf_ON = vip_hci.preproc.cube_px_resampling(psf_ON, psc_simu/psc_inst)
    psf_OFF = vip_hci.preproc.frame_px_resampling(psf_OFF, psc_simu/psc_inst)
    # clip negative values due to resampling
    psf_ON = psf_ON.clip(min=0)
    psf_OFF = psf_OFF.clip(min=0)
    
    """ calculate detector integration time (DIT) """
    DIT = adi_cube_duration/ncube
    
    """ rescale PSFs to stellar flux """
    # magnitude 5 star flux [e-/s], from Roy (Jan 8, 2019)
    mag5_ADU_all = {'Lp' : 1.834e+09,
                    'Mp' : 5.204e+08,
                    'N1': 2.291e+08,
                    'N2': 2.398e+08}
    # rescale to star flux
    star_flux = DIT * mag5_ADU_all[band] * 10**(-0.4*(mag - 5))
    psf_OFF *= star_flux
    psf_ON *= star_flux
    
    """ add background and photon noise ~ N(0,1) * sqrt(psf)"""
    if add_bckg is True:
        # background flux [e-/s/pix], from Roy (Jan 8, 2019)
        bckg_ADU_all = {'Lp' : 2.754e+05,
                        'Mp' : 2.010e+06,
                        'N1': 1.059e+08,
                        'N2': 3.293e+08}
        # add the background * transmission
        bckg_flux = DIT * bckg_ADU_all[band]
        psf_ON += bckg_flux*trans
        # add photon noise
        noise = np.random.standard_normal(psf_ON.shape) * np.sqrt(psf_ON)
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
    # VIP post-processing algorithm
    algo = vip_hci.medsub.median_sub
    # psf after post-processing
    out, derot, psf_pp = algo(psf_ON, pa, full_output=True)
    # contrast curve after post-processing
    cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop, fwhm, psc_inst/1e3, \
            starphot, algo=algo, nbranch=1, sigma=5, debug=False, plot=False)
    
    """ saving to fits files """
    fits.writeto(os.path.join(path_output, 'psf_' + filename + '.fits'), psf_pp, overwrite=True)
    hdu = fits.PrimaryHDU(cc_pp)
    hdu.writeto(os.path.join(path_output, 'cc_' + filename + '.fits'), overwrite=True)
    
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
        plt.xlim(left=0)
        plt.ylim([1e-7, 1e-0])
        plt.show(block=False)
        plt.savefig(os.path.join(path_output, 'cc_' + filename + '.png'), dpi=300, transparent=True)
