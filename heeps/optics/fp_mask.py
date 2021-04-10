from .lens import lens
from .vortex_init import vortex_init
from .lyotmask_init import lyotmask_init
import os.path
import numpy as np
from heeps.util.img_processing import resize_img, pad_img
from astropy.io import fits
import proper

def fp_mask(wf, mode='RAVC', vc_zoffset=0, verbose=False, **conf):

    # case 1: vortex coronagraphs
    if mode in ['CVC', 'RAVC']:
        if verbose is True:
            print('   apply vortex phase mask')                        
        # load vortex calibration files: psf_num, vvc, perf_num
        conf = vortex_init(verbose=verbose, **conf)
        # propagate to vortex
        lens(wf, offset_after=vc_zoffset, **conf)
        # apply vortex
        scale_psf = wf._wfarr[0,0]/conf['psf_num'][0,0]
        wf_corr = (conf['psf_num']*conf['vvc'] - conf['perf_num'])*scale_psf
        wf._wfarr = wf._wfarr*conf['vvc'] - wf_corr
        # propagate to lyot stop
        lens(wf, offset_before=-vc_zoffset, **conf)
    
    # case 2: classical Lyot
    elif mode in ['CLC']:
        if verbose is True:
            print('   apply classical lyot mask')
        # load lyotmask amplitude file
        conf = lyotmask_init(verbose=verbose, **conf)
        # propagate to lyot mask
        lens(wf, **conf)        
        # apply lyot mask
        wf._wfarr.real *= conf['lyotmask']
        # propagate to lyot stop
        lens(wf, **conf)

    return wf