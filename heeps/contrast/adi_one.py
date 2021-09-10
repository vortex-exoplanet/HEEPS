from heeps.util.save2fits import save2fits
from heeps.util.img_processing import crop_cube
from .background import background
import vip_hci
import numpy as np
from astropy.io import fits
import os.path
import warnings

def adi_one(dir_output='output_files', band='L', mode='RAVC', add_bckg=False,
        pscale=5.47, dit=0.3, mag=5, lat=-24.59, dec=-5, rim=19,
        app_strehl=0.64, nscreens=None, ndet=None, tag=None, student_distrib=True,
        savepsf=False, savefits=False, starphot=1e11, verbose=False, **conf):
    
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
        add_bckg (bool)
            true means background flux and photon noise are added 
        pscale (float):
            pixel scale in mas/pix (e.g. METIS LM=5.47, NQ=6.79)
        dit (float):
            detector integration time in s
        mag (float):
            star magnitude at selected band
        lat (float):
            telescope latitude in deg (Armazones=-24.59 ,Paranal -24.63)
        dec (float):
            star declination in deg (e.g. 51 Eri -2.47)
        rim (int):
            psf image radius in pixels
        app_strehl (float):
            APP Strehl ratio
        nscreens (int):
            number of screens of PSF cube taken into account, default to None for full cube
        ndet (int):
            size of the screens at the detector, default to None for full size
        tag (str):
            tag added to the saved filename
        student_distrib (bool):
            true if using a Student's t-distribution sensitivity, else Gaussian
        savepsf (bool):
            true if ADI psf is saved in a fits file
        savefits (bool):
            true if ADI contrast curve is saved in a fits file
        starphot (float):
            normalization factor for aperture photometry with VIP

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
    assert psf_ON.ndim == 3, "on-axis PSF cube must be 3-dimensional"
    # cut/crop cube
    if nscreens != None:
        psf_ON = psf_ON[:nscreens]
    if ndet != None:
        psf_ON = crop_cube(psf_ON, ndet)
    # get normalized off-axis PSF (planet)
    psf_OFF = fits.getdata(loadname%'offaxis')
    if psf_OFF.ndim == 3:
        psf_OFF = psf_OFF[0,:,:] # only first frame
    if verbose is True:
        print('Load PSFs for ADI')
        print('   mode=%s, band=%s'%(mode, band))
        print('   ncube=%s, ndet=%s'%(psf_ON.shape[0], psf_ON.shape[1]))
        print('   pscale=%s mas, dit=%s s'%(pscale, dit))
    # add background and photon noise: include star flux and HCI mode transmission
    if add_bckg is True:
        conf.update(mode=mode, dit=dit, mag=mag)
        psf_ON, psf_OFF = background(psf_ON, psf_OFF, verbose=True, **conf)
    # apply APP Strehl
    if 'APP' in mode:
        psf_OFF *= app_strehl

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
    ap_flux = vip_hci.metrics.aperture_flux(psf_OFF_crop, [rim], [rim], \
            fwhm, verbose=False)[0]
    # normalize to starphot (for VIP)
    if starphot is None:
        starphot = ap_flux
    else:
        psf_ON *= starphot/ap_flux
        psf_OFF_crop *= starphot/ap_flux

    """ parallactic angles for ADI """
    # hour angle to deg conversion
    ha = 360/24
    # angles in rad
    hr = np.deg2rad(np.linspace(-ha/2, ha/2, psf_ON.shape[0]))
    dr = np.deg2rad(dec)
    lr = np.deg2rad(lat)
    # parallactic angle in deg
    pa = -np.rad2deg(np.arctan2(-np.sin(hr), np.cos(dr)*np.tan(lr) \
            - np.sin(dr)*np.cos(hr)))
    
    """ VIP: post-processing (ADI, ADI-PCA,...) """
    # VIP post-processing algorithm
    algo = vip_hci.medsub.median_sub
    # contrast curve after post-processing (pscale in arcsec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # for AstropyDeprecationWarning
        cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop, \
                fwhm, pscale/1e3, starphot, algo=algo, nbranch=1, sigma=5, \
                debug=False, plot=False, verbose=verbose)
    # angular separations (in arcsec)
    sep = cc_pp.loc[:,'distance_arcsec'].values
    # sensitivities (Student's or Gaussian distribution)
    distrib = 'sensitivity_student' if student_distrib == True else 'sensitivity_gaussian'
    sen = cc_pp.loc[:,distrib].values
    # filename for fitsfiles
    if add_bckg is True:
        name = 'adi_bckg%s_mag%s'%(int(add_bckg), mag)
    else:
        name = 'adi_bckg%s'%int(add_bckg)
    # tag
    tag = '_%s'%tag.replace('/', '_') if tag != None else ''
    # save contrast curves as fits file
    if savefits == True:
        save2fits(np.array([sep,sen]), 'cc_%s%s%s'%(name, '_%s_%s', tag), dir_output=dir_output, band=band, mode=mode)
    # psf after post-processing
    if savepsf is True:
        _, _, psf_pp = algo(psf_ON, pa, full_output=True, verbose=False)
        save2fits(psf_pp, 'psf_%s%s%s'%(name, '_%s_%s', tag), dir_output=dir_output, band=band, mode=mode)

    return sep, sen
    

