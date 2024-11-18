import matplotlib.pyplot as plt
import proper
import numpy as np
from astropy.io import fits
import warnings

def get_wf(wf, plane, npupil=None, margin=0, savefits=False):

    assert 'amp' in plane or 'phi' in plane, 'plane must contain "amp" or "phi"'

    if 'amp' in plane:
        img = proper.prop_get_amplitude(wf)
    elif 'phi' in plane:
        img = proper.prop_get_phase(wf)
    ngrid = img.shape[0]
    if npupil is None:
        npupil = ngrid
    start = int((ngrid - npupil + npupil%2)/2 - margin)
    start = 0 if start < 0 else start
    end = ngrid - start + npupil%2
    img = img[start:end,start:end]
    if savefits is True:
        fits.writeto('%s_npupil=%s_margin=%s.fits'%(plane, npupil, margin),
            np.float32(img), overwrite=True)
    
    return img


def show_wf(wf, plane, npupil=None, margin=0, log=False, savefits=False):
    
    img = get_wf(wf, plane, npupil=npupil, margin=margin, savefits=savefits)
    if log is True:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") # divide by zero encountered in log10
            plt.figure(); plt.imshow(np.log10(img), origin='lower'); plt.colorbar();
    else:
        plt.figure(); plt.imshow(img, origin='lower'); plt.colorbar();

    return img