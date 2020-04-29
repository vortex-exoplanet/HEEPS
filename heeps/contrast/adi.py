import matplotlib.pyplot as plt
import vip_hci
import numpy as np
from astropy.io import fits
import os.path


def adi(path_offaxis='output_files', path_onaxis='output_files', 
        path_output='output_files', prefix='', scao_name='compass', 
        cube_duration=3600, cube_samp=300, adi_cube_duration=3600, 
        adi_cube_samp=0, adi_cube_avg=0, lat=-24.59, dec=-5, band='L', 
        mag=5, mode='CVC', pscale=5.47, rim=19, add_bckg=True, 
        tagname='', calc_trans=False, plot_cc=False, verbose=True):
    """ 
    This function calculates and draws the contrast curve (5-sigma sensitivity) 
    for a specific set of off-axis PSF and on-axis PSF (or cube of PSFs),
    using the VIP package to perform ADI post-processing.
    
    Args:
        path_offaxis (str):
            Path to off-axis PSFs
        path_onaxis (str):
            Path to on-axis (cubes of) PSFs
        path_output (str):
            Path to output files
        prefix (str):
            Optional loadfiles prefix
        scao_name (str):
            SCAO simulator name
        cube_duration (int):
            SCAO cube duration in seconds
        cube_samp (int):
            SCAO cube sampling in ms
        adi_cube_duration (int):
            ADI cube duration in seconds, default to 3600s (1h).
        adi_cube_samp (int):
            ADI cube sampling in ms. If 0 (default), same as cube_samp.
        adi_cube_avg (int):
            ADI cube averaging in ms. If 0 (default), no averaging.
        lat (float):
            Telescope latitude in deg (Armazones=-24.59 ,Paranal -24.63)
        dec (float):
            Star declination in deg (e.g. 51 Eri -2.47)
        band (str):
            Spectral band (e.g. 'L', 'M', 'N1', 'N2')
        mag (float):
            Star magnitude at band
        mode (str):
            HCI mode: ELT, VC, RAVC, APP, CL4, CL5,...
        pscale (float):
            Pixel scale in mas/pix (e.g. METIS LM=5.47, NQ=6.79)
        rim (int):
            Psf image radius in pixels
        add_bckg (bool)
            True means background flux and photon noise are added 
        calc_trans (bool):
            True means transmission is calculated from fits files
        plot_cc (bool):
            True means contrast curve is plotted
        verbose (bool):
            True means VPI functions print to stdout intermediate info and timing
        tagname (str):
            Tag name appended to savename, e.g. '_test'
    """
    # format input folders
    path_offaxis = os.path.normpath(os.path.expandvars(path_offaxis))
    path_onaxis = os.path.normpath(os.path.expandvars(path_onaxis))
    path_output = os.path.normpath(os.path.expandvars(path_output))
    
    """ filenames """
    loadname = '%s%s_%s_%s_%s.fits'%(prefix,'%s','%s',band,'%s')
    savename = '%s%ss_samp%sms_dec%s_%smag%s_bckg%s_%s%s' \
            %(scao_name, cube_duration, cube_samp, dec, band, mag, int(add_bckg), mode, tagname)
    
    """ transmission : ratio of intensities (squared amplitudes) in Lyot-Stop plane """
    if calc_trans is True:
        I_ELT = fits.getdata(os.path.join(path_offaxis, \
                loadname%('offaxis', 'LS', 'ELT')))**2
        I_OFFAXIS = fits.getdata(os.path.join(path_offaxis, \
                loadname%('offaxis', 'LS', mode)))**2
        trans = np.sum(I_OFFAXIS)/np.sum(I_ELT)
    else:
        trans_all = {'ELT': 1.,
                    'RAVC': 0.3392759549914341,
                     'CVC': 0.9012406091763115,
                     'APP': 1.,
                     'CLC': 0.502257180317745} # CLC occulter diam = 4lam/D
        trans = trans_all[mode]
    
    """ get normalized off-axis PSF (single) """
    # total flux of the non-coronagraphic PSF
    psf_ELT = fits.getdata(os.path.join(path_offaxis, \
                loadname%('offaxis', 'PSF', 'ELT')))
    if psf_ELT.ndim == 3: # if cube
        psf_ELT = psf_ELT[0,:,:]
    ELT_flux = np.sum(psf_ELT)
    # normalized off-axis PSF
    psf_OFF = fits.getdata(os.path.join(path_offaxis, \
                loadname%('offaxis', 'PSF', mode)))
    if psf_OFF.ndim == 3: # if cube
        psf_OFF = psf_OFF[0,:,:]
    psf_OFF /= ELT_flux
    
    """ get normalized on-axis PSFs (cube, resampled, and averaged) """
    # load cube, and format to 3D
    psf_ON = fits.getdata(os.path.join(path_onaxis, \
                loadname%('onaxis', 'PSF', mode)))
    if psf_ON.ndim != 3:
        psf_ON = np.array(psf_ON, ndmin=3)
    # save PSF initial shape
    (ncube, xon, yon) = psf_ON.shape
    # resample based on ADI sampling vs simulation sampling
    if adi_cube_samp == 0:
        adi_cube_samp = cube_samp
    nsamp = max(1, int(adi_cube_samp/cube_samp + .5))
    psf_ON = psf_ON[::nsamp,:]
    # average frames
    navg = max(1, int(adi_cube_avg/adi_cube_samp + .5))
    if navg > 1:
        end = ncube - ncube % navg
        psf_ON = np.mean(psf_ON[:end].reshape(-1, navg, xon, yon), axis=1)
    # normalized coronagraphic (on-axis) PSFs
    psf_ON /= ELT_flux
    ncube = psf_ON.shape[0]
    
    """ calculate detector integration time (DIT) """
    DIT = adi_cube_duration/ncube
    
    """ rescale PSFs to stellar flux """
    # magnitude 0 star flux [e-/s], from Roy (Jan 21, 2020)
    mag_ref = 0
    star_flux_all = {'L' : 8.999e+10, #HCI-L long
                     'M' : 2.452e+10, #CO ref
                     'N1': 3.684e+10, #GeoSnap N1
                     'N2': 3.695e+10, #GeoSnap N2
                    'N1a': 2.979e+10, #Aquarius N1
                    'N2a': 2.823e+10} #Aquarius N2
    # rescale to star flux
    star_flux = DIT * star_flux_all[band] * 10**(-0.4*(mag - mag_ref))
    if mode == 'APP':
        star_flux *= 0.48 # for APP separated PSFs (4% leakage term)
    psf_OFF *= star_flux
    psf_ON *= star_flux
    if mode == 'APP':
        psf_OFF *= 0.64 # for APP 64% Strehl ratio
    
    """ add background and photon noise ~ N(0,1) * sqrt(psf) """
    if add_bckg is True:
        # background flux [e-/s/pix], from Roy (Jan 21, 2020)
        bckg_flux_all = {'L' : 8.878e+04, #HCI-L long
                         'M' : 6.714e+05, #CO ref
                         'N1': 4.725e+07, #GeoSnap N1
                         'N2': 1.122e+08, #GeoSnap N2
                        'N1a': 9.630e+07, #Aquarius N1
                        'N2a': 2.142e+08} #Aquarius N2
        # add the background * transmission
        bckg_flux = DIT * bckg_flux_all[band]
        psf_ON += bckg_flux*trans
        # add photon noise
        noise = np.random.standard_normal(psf_ON.shape) * np.sqrt(psf_ON)
        psf_ON += noise
    
    """ VIP: aperture photometry of psf_OFF used to scale the contrast """
    # get the center pixel
    (xoff, yoff) = psf_OFF.shape
    (cx, cy) = (int(xoff/2), int(yoff/2))
    # fit a 2D Gaussian --> output: fwhm, x-y centroid
    fit = vip_hci.var.fit_2dgaussian(psf_OFF[cx-rim:cx+rim+1, \
            cy-rim:cy+rim+1], True, (rim,rim), debug=False, full_output=True)
    # derive the FWHM
    fwhm = np.mean([fit['fwhm_x'],fit['fwhm_y']])
    # recenter and crop
    shiftx, shifty = rim-fit['centroid_x'], rim-fit['centroid_y']
    psf_OFF = vip_hci.preproc.frame_shift(psf_OFF, shiftx, shifty)
    psf_OFF_crop = psf_OFF[cx-rim:cx+rim+1, cy-rim:cy+rim+1]
    # FWHM aperture photometry of psf_OFF_crop
    starphot = vip_hci.metrics.aperture_flux(psf_OFF_crop, [rim], [rim], \
            fwhm, verbose=False)[0]
    
    """ parallactic angles for ADI """
    # duration -> hour angle conversion
    ha = adi_cube_duration/3600/24*360
    # angles in rad
    hr = np.deg2rad(np.linspace(-ha/2, ha/2, ncube))
    dr = np.deg2rad(dec)
    lr = np.deg2rad(lat)
    # parallactic angle in deg
    pa = -np.rad2deg(np.arctan2(-np.sin(hr), np.cos(dr)*np.tan(lr) \
            - np.sin(dr)*np.cos(hr)))
    
    """ VIP: post-processing (ADI, ADI-PCA,...) """
    # VIP post-processing algorithm
    algo = vip_hci.medsub.median_sub
    # psf after post-processing
    if False:
        out, derot, psf_pp = algo(psf_ON, pa, full_output=True, verbose=False)
        fits.writeto(os.path.join(path_output, 'psf_' + savename + '.fits'), \
                psf_pp, overwrite=True)
    # contrast curve after post-processing (pscale in arcsec)
    cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop, \
            fwhm, pscale/1e3, starphot, algo=algo, nbranch=1, sigma=5, \
            debug=False, plot=False, verbose=verbose)
    hdu = fits.PrimaryHDU(cc_pp)
    hdu.writeto(os.path.join(path_output, 'cc_' + savename + '.fits'), overwrite=True)
    
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
        plt.savefig(os.path.join(path_output, 'cc_' + savename + '.png'), \
                dpi=300, transparent=True)
