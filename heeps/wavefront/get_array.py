import matplotlib.pyplot as plt
import proper
import numpy as np
from astropy.io import fits

def get_wf(wf, plane, npupil=None, margin=0, savefits=False):

    assert plane in ['amp','phi'], 'plane must be "amp" or "phi"'

    if plane is 'amp':
        img = proper.prop_get_amplitude(wf)
    elif plane is 'phi':
        img = proper.prop_get_phase(wf)
    ngrid = img.shape[0]
    if npupil is None:
        npupil = ngrid
    start = int((ngrid - npupil + npupil%2)/2 - margin)
    start = 0 if start < 0 else start
    end = ngrid - start + npupil%2
    img = img[start:end,start:end]
    if savefits is True:
        fits.writeto('%s_npupil=%s_margin=%s.fits'%(plane, npupil, margin), \
            np.float32(img), overwrite=True)
    
    return img


def show_wf(wf, plane, npupil=None, margin=0, savefits=False):
    
    img = get_wf(wf, plane, npupil=npupil, margin=margin, savefits=savefits)
    plt.imshow(img, origin='lower')

    return img