import numpy as np
from astropy.io import fits

def background(psf_ON, psf_OFF, mode='RAVC', lam=3.8e-6, cube_duration=3600, 
        mag=5, mag_ref=0, flux_star=9e10, flux_bckg=9e4, app_single_psf=0.48, 
        f_vc_trans=None, f_app_trans=None, seed=123456, verbose=False, **conf):

    """ 
    This function applies background and photon noise to intup PSFs (off-axis and on-axis), 
    incuding transmission, star flux, and components transmittance.

    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        mode (str):
            HCI mode: RAVC, CVC, APP, CLC            
        lam (float):
            wazvelength in m
        cube_duration (int):
            ADI sequence duration in s
        mag (float):
            star magnitude
        mag_ref (float):
            reference magnitude for star and background fluxes
        flux_star (float):
            star flux at reference magnitude
        flux_bckg (float):
            background flux at reference magnitude
        app_single_psf (float):
            APP single PSF (4% leakage)
        f_vc_trans (str):
            path to VC transmittance fits file
        f_app_trans
            path to APP transmittance fits file
        seed (int):
            seed used by numpy.random process
    
    Return:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
    """

    # calculate offaxis-PSF transmission
    offaxis_trans = np.sum(psf_OFF)
    # apply mask transmittance
    data = None
    if 'VC' in mode:
        data = fits.getdata(f_vc_trans)
    elif 'APP' in mode:
        data = fits.getdata(f_app_trans)
    if data is not None:
        mask_trans = np.interp(lam*1e6, data[0], data[1])
        psf_OFF *= mask_trans
        psf_ON *= mask_trans
    # apply correction for APP single PSF
    if 'APP' in mode:
        psf_OFF *= app_single_psf
        psf_ON *= app_single_psf
    # detector integration time
    DIT = cube_duration/psf_ON.shape[0]
    # rescale PSFs to star signal
    star_signal = DIT * flux_star * 10**(-0.4*(mag - mag_ref))
    psf_OFF *= star_signal
    psf_ON *= star_signal
    # add background
    bckg_noise = DIT * flux_bckg * offaxis_trans * mask_trans
    psf_ON += bckg_noise
    # add photon noise ~ N(0, sqrt(psf))
    np.random.seed(seed)
    psf_ON += np.random.normal(0, np.sqrt(psf_ON))

    if verbose is True:
        print('   offaxis_trans=%3.4f, mask_trans=%3.4f,'%(offaxis_trans, mask_trans))
        print('   mag=%s, DIT=%3.3f'%(mag, DIT))
        print('   star_signal=%3.2E, bckg_noise=%3.2E'%(star_signal, bckg_noise))

    return psf_ON, psf_OFF