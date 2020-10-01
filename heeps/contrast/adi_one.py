from heeps.util.save2fits import save2fits
import vip_hci
import numpy as np
from astropy.io import fits
import os.path

def adi_one(dir_output='output_files', band='L', mode='RAVC', mag=5, mag_ref=0, 
        flux_star=9e10, flux_bckg=9e4, add_bckg=True, pscale=5.47, 
        cube_duration=3600, lat=-24.59, dec=-5, rim=19, app_strehl=0.64, 
        app_single_psf=0.48, student_distrib=True, seed=123456, savepsf=False, 
        savefits=False, verbose=False, **conf):
    
    """ 
    This function calculates and draws the contrast curve (5-sigma sensitivity) 
    for a specific set of off-axis PSF and on-axis PSF (or cube of PSFs),
    using the VIP package to perform ADI post-processing.
    
    Args:
        dir_output (str):
            path to output fits files and input PSFs (onaxis & offaxis)
        band (str):
            spectral band (e.g. 'L', 'M', 'N1', 'N2')
        mode (str):
            HCI mode: RAVC, CVC, APP, CLC
        mag (float):
            star magnitude at selected band
        mag_ref (float):
            reference magnitude for star and background fluxes
        flux_star (float):
            star flux at reference magnitude
        flux_bckg (float):
            background flux at reference magnitude
        add_bckg (bool)
            true means background flux and photon noise are added 
        pscale (float):
            pixel scale in mas/pix (e.g. METIS LM=5.47, NQ=6.79)
        cube_duration (int):
            cube duration in seconds, default to 3600s (1h).
        lat (float):
            telescope latitude in deg (Armazones=-24.59 ,Paranal -24.63)
        dec (float):
            star declination in deg (e.g. 51 Eri -2.47)
        rim (int):
            psf image radius in pixels
        app_strehl (float):
            APP Strehl ratio
        app_single_psf (float):
            APP single PSF (4% leakage)
        student_distrib (bool):
            true if using a Student's t-distribution sensitivity, else Gaussian

    Return:
        sep (float ndarray):
            angular separation in arcsec
        sen (float ndarray):
            5-sigma sensitivity (contrast)

    """

    # PSF filenames
    loadname = os.path.join(dir_output, '%s_PSF_%s_%s.fits'%('%s', band, mode))
    # get normalized on-axis PSFs (star)
    psf_ON = fits.getdata(loadname%'onaxis')
    if psf_ON.ndim != 3: # must be a cube (3D)
        psf_ON = np.array(psf_ON, ndmin=3)
    ncube = psf_ON.shape[0]
    # get normalized off-axis PSF (planet)
    psf_OFF = fits.getdata(loadname%'offaxis')
    if psf_OFF.ndim == 3: # only one frame
        psf_OFF = psf_OFF[0,:,:]
    # calculate transmission
    trans = np.sum(psf_OFF)
    # apply correction for APP PSFs
    if mode == 'APP':
        psf_OFF *= app_single_psf*app_strehl
        psf_ON *= app_single_psf
    # detector integration time (DIT)
    DIT = cube_duration/ncube
    # rescale PSFs to star signal
    star_signal = DIT * flux_star * 10**(-0.4*(mag - mag_ref))
    psf_OFF *= star_signal
    psf_ON *= star_signal
    if verbose is True:
        print('Load PSFs for ADI')
        print('   mode=%s, band=%s, pscale=%s'%(mode, band, pscale))
        print('   DIT=%3.3f, mag=%s, star_signal=%3.2E'%(DIT, mag, star_signal))
    # add background and photon noise ~ N(0, sqrt(psf))
    if add_bckg is True:            
        bckg_noise = DIT * flux_bckg * trans
        psf_ON += bckg_noise
        np.random.seed(seed)
        phot_noise = np.random.normal(0, np.sqrt(psf_ON))
        psf_ON += phot_noise
        if verbose is True:
            print('   add_bckg=%s, trans=%3.4f, bckg_noise=%3.2E'\
                %(add_bckg, trans, bckg_noise))

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
    ha = cube_duration/3600/24*360
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
    # contrast curve after post-processing (pscale in arcsec)
    cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop, \
            fwhm, pscale/1e3, starphot, algo=algo, nbranch=1, sigma=5, \
            debug=False, plot=False, verbose=verbose)
    # angular separations (in arcsec)
    sep = cc_pp.loc[:,'distance_arcsec'].values
    # sensitivities (Student's or Gaussian distribution)
    distrib = 'sensitivity_student' if student_distrib == True else 'sensitivity_gaussian'
    sen = cc_pp.loc[:,distrib].values

    # save contrast curves as fits file
    if savefits == True:
        name = 'cc_adi_bckg%s_mag%s'%(int(add_bckg), mag)
        save2fits(np.array([sep,sen]), name, dir_output=dir_output, band=band, mode=mode)

    # psf after post-processing
    if savepsf is True:
        out, derot, psf_pp = algo(psf_ON, pa, full_output=True, verbose=False)
        name = 'psf_adi_bckg%s_mag%s'%(int(add_bckg), mag)
        save2fits(psf_pp, name, dir_output=dir_output, band=band, mode=mode)

    return sep, sen
    

