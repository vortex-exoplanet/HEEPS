from .background import background
from heeps.util.save2fits import save2fits
from heeps.util.img_processing import crop_cube
from heeps.util.paralang import paralang
from heeps.util.psf_template import psf_template
import multiprocessing as mpro
import vip_hci
import numpy as np
from astropy.io import fits
import os.path
import warnings
from packaging import version

vvip = vip_hci.__version__

def cc_adi(dir_output='output_files', band='L', mode='RAVC', add_bckg=False,
        pscale=5.47, dit=0.3, mag=5, lat=-24.59, dec=-5, app_strehl=0.64, 
        nscreens=None, ndet=None, tag=None, f_oat=None, student_distrib=True, 
        savepsf=False, starphot=1e11, duration=3600, cpu_count=None,
        savefits=False, verbose=False, **conf):

    """
    This function calculates and draws the contrast curve (5-sigma sensitivity) 
    for a specific set of off-axis PSF and on-axis cube of PSFs,
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
        app_strehl (float):
            APP Strehl ratio
        nscreens (int):
            number of screens of PSF cube taken into account, default to None for full cube
        ndet (int):
            size of the screens at the detector, default to None for full size
        tag (str):
            tag added to the saved filename
        f_oat (str):
            path to off-axis transmission file, default to None
        student_distrib (bool):
            true if using a Student's t-distribution sensitivity, else Gaussian
        savepsf (bool):
            true if ADI psf is saved in a fits file
        savefits (bool):
            true if ADI contrast curve is saved in a fits file
        starphot (float):
            normalization factor for aperture photometry with VIP
        cpu_count (int):
            number of CPU cores to use, 'None' means use maximum number of cores


    Return:
        sep (float ndarray):
            angular separation in arcsec
        adi (float ndarray):
            ADI 5-sigma sensitivity (contrast)
    """

    # load PSFs: on-axis (star) and off-axis (planet)
    loadname = os.path.join(dir_output, '%s_PSF_%s_%s.fits'%('%s', band, mode))
    psf_OFF = fits.getdata(loadname%'offaxis')
    psf_ON = fits.getdata(loadname%'onaxis')
    header_ON = fits.getheader(loadname%'onaxis')
    assert psf_ON.ndim == 3, "on-axis PSF cube must be 3-dimensional"
    assert psf_OFF.ndim == 2, "off-axis PSF frame must be 2-dimensional"
    # cut/crop cube
    if nscreens is not None:
        psf_ON = psf_ON[:nscreens]
    if ndet is not None:
        psf_ON = crop_cube(psf_ON, ndet)
    if verbose is True:
        print('Apply ADI technique: add_bckg=%s'%add_bckg)
        print('\u203e'*20)
        print('   mode=%s, band=%s, pscale=%s mas'%(mode, band, pscale))
        print('   ncube=%s, ndet=%s, duration=%s s'%(psf_ON.shape[0], psf_ON.shape[1], duration))
    # add background and photon noise: include star flux and HCI mode transmission
    if add_bckg is True:
        conf.update(mode=mode, dit=dit, mag=mag)
        psf_ON, psf_OFF = background(psf_ON, psf_OFF, header=header_ON,
            verbose=True, **conf)
    # apply APP Strehl
    if 'APP' in mode:
        psf_OFF *= app_strehl
    # parallactic angle in deg
    pa = paralang(psf_ON.shape[0], dec, lat, duration=duration)
    # get off-axis transmission
    if 'VC' in mode and f_oat is not None:
        OAT = fits.getdata(f_oat)
        if version.parse(vvip) < version.parse("1.3.0"):
            OAT = (OAT[1], OAT[0])
        else:
            pass
        if verbose is True:
            print("   load vortex off-axis transmission from '%s'"%os.path.basename(f_oat))
    else:
        OAT = None
    # aperture photometry of an off-axis PSF template, used to scale the contrast
    psf_OFF_crop, fwhm, ap_flux = psf_template(psf_OFF)
    # normalize to starphot (for VIP)
    if starphot is None:
        starphot = ap_flux
    else:
        psf_ON *= starphot/ap_flux
        psf_OFF_crop *= starphot/ap_flux
    # VIP post-processing algorithm
    if cpu_count == None:
        cpu_count = mpro.cpu_count()
    algo_dict = dict(nproc=cpu_count)
    if version.parse(vvip) <= version.parse("1.0.3"):
        algo = vip_hci.medsub.median_sub
    else:
        algo = vip_hci.psfsub.median_sub    
    # contrast curve after post-processing (pscale in arcsec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # for AstropyDeprecationWarning
        cc_pp = vip_hci.metrics.contrast_curve(psf_ON, pa, psf_OFF_crop,
                fwhm, pscale/1e3, starphot, algo=algo, nbranch=1, sigma=5,
                debug=False, plot=False, transmission=OAT, imlib='opencv',
                verbose=verbose, **algo_dict)
        if verbose is True:
            print('')
    # angular separations (in arcsec)
    sep = cc_pp.loc[:,'distance_arcsec'].values
    # sensitivities (Student's or Gaussian distribution)
    distrib = 'sensitivity_student' if student_distrib == True else 'sensitivity_gaussian'
    adi = cc_pp.loc[:,distrib].values
    # filename for fitsfiles
    if add_bckg is True:
        name = 'adi_bckg%s_mag%s'%(int(add_bckg), mag)
    else:
        name = 'adi_bckg%s'%int(add_bckg)
    # tag
    tag = '_%s'%tag.replace('/', '_') if tag != None else ''
    # save contrast curves as fits file
    if savefits == True:
        save2fits(np.array([sep, adi]), 'cc_%s%s%s'%(name, '_%s_%s', tag),
            dir_output=dir_output, band=band, mode=mode)
    # psf after post-processing
    if savepsf is True:
        _, _, psf_pp = algo(psf_ON, pa, full_output=True, verbose=False)
        save2fits(psf_pp, 'psf_%s%s%s'%(name, '_%s_%s', tag), 
            dir_output=dir_output, band=band, mode=mode)

    return sep, adi