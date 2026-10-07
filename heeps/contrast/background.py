import numpy as np
from astropy.io import fits

def background(psf_ON, psf_OFF, header=None, mode='RAVC', lam=3.8e-6, dit=0.3, 
        mag=5, mag_ref=0, flux_star=9e10, flux_bckg=9e4, app_single_psf=0.48, 
        f_vc_trans=None, f_app_trans=None, seed=123456, verbose=False, 
        call_ScopeSim=False, **conf):

    """
    This function applies background and photon noise to intup PSFs (off-axis
    and on-axis), including transmission, star flux, and components transmittance.

    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        mode (str):
            HCI mode: RAVC, CVC, APP, CLC
        lam (float):
            wavelength in m
        dit (float):
            detector integration time in s
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
        call_ScopeSim (bool):
            true if interfacing ScopeSim

    Return:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
    """

    # calculate offaxis-PSF transmission
    thruput = np.sum(psf_OFF)
    # load mask transmittance
    if 'VC' in mode:
        data = fits.getdata(f_vc_trans)
        mask_trans = np.interp(lam*1e6, data[0], data[1])
    elif 'APP' in mode:
        data = fits.getdata(f_app_trans)
        mask_trans = np.interp(lam*1e6, data[0], data[1])
    else:
        mask_trans = 1
    # apply correction for APP single PSF (~48%)
    if 'APP' in mode:
        psf_OFF *= app_single_psf
        psf_ON *= app_single_psf

    # scopesim-heeps interface
    if call_ScopeSim:
        from heeps.contrast.sim_heeps import sim_heeps
        psf_ON, psf_OFF = sim_heeps(psf_ON, psf_OFF, header, **conf)
    else:
        # rescale PSFs to star signal
        star_signal = dit * flux_star * 10**(-0.4*(mag - mag_ref))
        psf_OFF *= star_signal * mask_trans
        psf_ON *= star_signal * mask_trans
        # add background
        bckg_noise = dit * flux_bckg * thruput * mask_trans
        psf_ON += bckg_noise
        # add photon noise ~ N(0, sqrt(psf))
        # default_rng with standard_normal is much faster than
        # np.random.normal(0, sigma_array) for large cubes, because
        # standard_normal draws from N(0,1) using a fast vectorised path
        # and then scale by sqrt(psf_ON) in a single multiply pass,
        # avoiding per-element sigma sampling
        rng = np.random.default_rng(seed)
        psf_sqrt = np.sqrt(psf_ON)
        psf_sqrt *= rng.standard_normal(psf_ON.shape)
        psf_ON += psf_sqrt

    if verbose:
        print('   dit=%s s, thruput=%.4f, mask_trans=%.4f,'%(dit, thruput, mask_trans))
        print('   mag=%s, star_signal=%.2e, bckg_noise=%.2e'%(mag, star_signal, bckg_noise))

    return psf_ON, psf_OFF