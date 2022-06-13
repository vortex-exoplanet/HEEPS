import numpy as np
import vip_hci

def psf_off_template(psf_OFF, recenter=True, ncrop=4):

    # get the center pixel (PSF must be centered in the image)
    (yoff, xoff) = psf_OFF.shape
    (cx, cy) = (xoff//2, yoff//2)
    # fit a 2D Gaussian --> output: fwhm, x-y centroid
    fit = vip_hci.var.fit_2dgaussian(psf_OFF, True, (cx, cy), debug=False, full_output=True)
    # derive the FWHM
    fwhm = np.mean([fit['fwhm_x'], fit['fwhm_y']])
    # recenter
    if recenter is True:
        yerr, xerr = fit['centroid_y_err'].values[0], fit['centroid_x_err'].values[0]
        assert (yerr < 0.5) and (xerr < 0.5), 'centroid (x,y) error = '\
                                    '(%.1e, %.1e), must be < 0.5 pixel'%(xerr, yerr)
        shiftx, shifty = cx-fit['centroid_x'].values[0], cy-fit['centroid_y'].values[0]
        psf_OFF = vip_hci.preproc.frame_shift(psf_OFF, shifty, shiftx)
    # FWHM aperture photometry
    ap_flux = vip_hci.metrics.aperture_flux(psf_OFF, [cy], [cx], fwhm, verbose=False)[0]
    # image radius (cropped), default to 4x FWHM
    rim = round(ncrop*fwhm)
    psf_OFF_crop = psf_OFF[cy-rim:cy+rim+1, cx-rim:cx+rim+1]

    return psf_OFF_crop, fwhm, ap_flux