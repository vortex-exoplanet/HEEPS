import numpy as np

from vip_hci.metrics import aperture_flux
from vip_hci.preproc import frame_shift
from vip_hci.var import fit_2dgaussian


def psf_template(psf, center=None, recenter=True, ncrop=4):

    # get the center pixel (PSF must be centered in the image)
    if center is None:
        (ny, nx) = psf.shape
        (cx, cy) = (nx//2, ny//2)
    else:
        (cx, cy) = center
    # fit a 2D Gaussian --> output: fwhm, x-y centroid
    fit = fit_2dgaussian(psf, True, (cx, cy), debug=False, full_output=True)
    # derive the FWHM
    fwhm = np.mean([fit['fwhm_x'], fit['fwhm_y']])
    # recenter
    if recenter is True:
        yerr, xerr = fit['centroid_y_err'].values[0], fit['centroid_x_err'].values[0]
        assert (yerr < 0.5) and (xerr < 0.5), 'centroid (x,y) error = '\
                                    '(%.1e, %.1e), must be < 0.5 pixel'%(xerr, yerr)
        shiftx, shifty = cx-fit['centroid_x'].values[0], cy-fit['centroid_y'].values[0]
        psf = frame_shift(psf, shifty, shiftx)
    # FWHM aperture photometry
    ap_flux = aperture_flux(psf, [cy], [cx], fwhm, verbose=False)[0]
    # image radius (cropped), default to 4x FWHM
    rim = round(ncrop*fwhm)
    psf_crop = psf[cy-rim:cy+rim+1, cx-rim:cx+rim+1]

    return psf_crop, fwhm, ap_flux
