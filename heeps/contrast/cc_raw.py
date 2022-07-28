from heeps.util.img_processing import crop_cube, get_radial_profile
from heeps.util.save2fits import save2fits
import os.path
from astropy.io import fits
import numpy as np


def cc_raw(dir_output='output_files', band='L', mode='RAVC', pscale=5.47,
        nscreens=None, ndet=None,tag=None, savefits=False, verbose=False, **conf):
    """
    """

    # load PSFs: on-axis (star) and off-axis (planet)
    loadname = os.path.join(dir_output, '%s_PSF_%s_%s.fits'%('%s', band, mode))
    psf_OFF = fits.getdata(loadname%'offaxis')
    psf_ON = fits.getdata(loadname%'onaxis')
    assert psf_ON.ndim in [2, 3], "on-axis PSF can be 2- or 3-dimensional"
    assert psf_OFF.ndim == 2, "off-axis PSF frame must be 2-dimensional"
    # force 3D
    psf_ON = np.array(psf_ON, ndmin=3)
    # cut/crop cube
    if nscreens is not None:
        psf_ON = psf_ON[:nscreens]
    if ndet is not None:
        psf_ON = crop_cube(psf_ON, ndet)
    if verbose is True:
        print('Raw contrast curve:')
        print('\u203e'*19)
        print('   mode=%s, band=%s, pscale=%s mas'%(mode, band, pscale))
        print('   ncube=%s, ndet=%s\n'%(psf_ON.shape[0], psf_ON.shape[1]))
    # average
    psf_ON_avg = np.mean(psf_ON, 0)
    # image radius
    rim = ndet // 2
    # angular separation in arcsec
    sep = pscale*1e-3*np.arange(rim)
    # radial profiles
    off = get_radial_profile(psf_OFF, (rim,rim), 1)[:-1]
    raw = get_radial_profile(psf_ON_avg, (rim,rim), 1)[:-1]
    # normalize by the peak of the off-axis PSF radial profile
    raw /= np.max(off)
    # tag
    tag = '_%s'%tag.replace('/', '_') if tag != None else ''
    # save contrast curves as fits file
    if savefits == True:
        save2fits(np.array([sep, raw]), 'cc_%s%s%s'%('raw', '_%s_%s', tag),
            dir_output=dir_output, band=band, mode=mode)

    return sep, raw