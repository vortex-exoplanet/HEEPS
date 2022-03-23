from heeps.optics.lens import lens
import proper
from astropy.io import fits
import os
import numpy as np

def detector(wf, ngrid=1024, ndet=365, dir_output='output_files', savefits=False, verbose=False, **conf):
    
    "Extract PSF on the detector."

    assert(ngrid >= ndet), 'Error: final image is bigger than initial grid size'

    # propagate to detector
    lens(wf, **conf)
    # get intensity (A^2)
    (psf, _) = proper.prop_end(wf, NOABS = False)
    # crop to detector size
    start = int(ngrid/2 - ndet/2) + 1
    end = int(ngrid/2 + ndet/2) + 1
    psf = psf[start:end, start:end]

    if verbose is True:
        print("   extract PSF on the detector: ndet=%s"%ndet)

    # save psf as fits file
    if savefits == True:
        os.makedirs(dir_output, exist_ok=True)
        filename = os.path.join(dir_output, 'PSF_IMG_%s.fits'%conf['band'])
        fits.writeto(filename, np.float32(psf), overwrite=True)

    return psf