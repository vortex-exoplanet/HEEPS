from photutils import CircularAperture, aperture_photometry
import astropy.units as u
import numpy as np


def get_psf_flux(psf, lam=3.81e-6, diam_ext=36.905, diam_nominal=38.542,
        ls_dRext=0.0209, pscale=5.47, verbose=False, **conf):

    # full width half max (= lam/D)
    fwhm = get_lamD_pix(lam=lam, diam_ext=diam_ext, diam_nominal=diam_nominal, 
        ls_dRextls_dRext=ls_dRext, pscale=pscale) 
    # photutils aperture photometry
    nimg = psf.shape[-1]
    aper = CircularAperture((nimg//2, nimg//2), r=fwhm/2)
    psf_flux = aperture_photometry(psf, aper)['aperture_sum'].data
    if verbose is True:
        print('photutils aperture photometry: psf_flux=%s'%np.round(psf_flux,5))

    return psf_flux


def get_lamD_pix(lam=3.81e-6, diam_ext=36.905, diam_nominal=38.542,
        ls_dRext=0.0209, pscale=5.47, **conf):
    
    # 1 lambda/D in mas
    lamD = get_lamD_mas(lam=lam, diam_ext=diam_ext, diam_nominal=diam_nominal, 
        ls_dRextls_dRext=ls_dRext)
    # 1 lambda/D in pixels
    lamD /= pscale
    
    return lamD


def get_lamD_mas(lam=3.81e-6, diam_ext=36.905, diam_nominal=38.542,
        ls_dRext=0.0209, **conf):
    
    # 1 lambda/D in mas
    diamLS = diam_ext - diam_nominal*ls_dRext
    lamD = lam/diamLS*u.rad.to('mas')

    return lamD