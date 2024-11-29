from .create_stop import create_stop
from heeps.util.img_processing import resize_img, pad_img
from heeps.util.save2fits import save2fits
import proper
import numpy as np
import os.path
from astropy.io import fits 


def pupil(pup=None, f_pupil='', lam=3.8e-6, ngrid=1024, npupil=285,
        pupil_img_size=40, diam_nominal=38, norm_I=True, savefits=False,
        verbose=False, **conf):
    """Creates a wavefront object at the entrance pupil plane. 
    The pupil is either loaded from a fits file, or created using 
    pupil parameters.

    Args:
        f_pupil: str
            path to a pupil fits file
        lam: float
            wavelength in m
        ngrid: int
            number of pixels of the wavefront array
        npupil: int
            number of pixels of the pupil
        pupil_img_size: float
            pupil image (for PROPER) in m
        diam_nominal: float
            nominal diameter in m
        norm_I: bool
            True if normalized intensity
    
    """

    if verbose is True:
        print("Entrance pupil:", end=' ')

    # initialize wavefront using PROPER
    beam_ratio = (diam_nominal/pupil_img_size) * (npupil/ngrid)
    wf = proper.prop_begin(diam_nominal, lam, ngrid, beam_ratio)

    # case 1: load pupil from data
    if pup is not None:
        if verbose is True:
            print("loaded from 'pup' array")
        pup = resize_img(pup, npupil)

    # case 2: load pupil from file
    elif os.path.isfile(f_pupil):
        if verbose is True:
            print("loaded from '%s'"%os.path.basename(f_pupil))
        pup = resize_img(fits.getdata(f_pupil), npupil)
    
    # case 3: create a pupil
    else:
        if verbose is True:
            print("created with", end=' ')
        conf.update(npupil=npupil,
                    pupil_img_size=pupil_img_size,
                    diam_nominal=diam_nominal)
        pup = create_stop(**conf)

    # normalize the entrance pupil intensity (total flux = 1)
    if norm_I is True:
        I_pup = pup**2
        pup = np.sqrt(I_pup/np.sum(I_pup))

    # save pupil as fits file
    if savefits == True:
        save2fits(pup, 'pupil', **conf)
    
    # pad with zeros and add to wavefront
    proper.prop_multiply(wf, pad_img(pup, ngrid))
    
    if verbose is True:
        print('\u203e'*15)

    return wf