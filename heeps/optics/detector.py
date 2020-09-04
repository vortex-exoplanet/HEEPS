import proper
from heeps.optics.lens import lens

def detector(wf, conf, verbose=False):
    f_lens = conf['focal']
    ndet = conf['ndet']
    ngrid = conf['ngrid']
    
    assert(ngrid >= ndet), 'Error: final image is bigger than initial grid size'

    # propagate to detector
    lens(wf, f_lens)
    # get intensity (A^2)
    (psf, _) = proper.prop_end(wf, NOABS = False)
    # crop to detector size
    start = int(ngrid/2 - ndet/2) + 1
    end = int(ngrid/2 + ndet/2) + 1
    psf = psf[start:end, start:end]

    if verbose is True:
        print("Extract PSF on the detector\n")

    return psf
