from .lens import lens
from .vortex_init import vortex_init
from .lyotmask_init import lyotmask_init
import numpy as np
import proper

def fp_mask(wf, mode='RAVC', vc_zoffset=0, vc_chrom_leak=2e-3,
        lt_dist_start=0, lt_diam_start=0, lt_dist_end=0, lt_diam_end=0,
        diam_vpm=16273, diam_lom=16273, add_cl_vort=False, verbose=False, **conf):

    # case 1: vortex coronagraphs
    if mode in ['CVC', 'RAVC']:
        if verbose is True:
            print('   apply vortex phase mask')                        
        # load vortex calibration files: psf_num, vvc, perf_num
        conf = vortex_init(verbose=verbose, **conf)
        # propagate to vortex
        lens(wf, zoffset_before=vc_zoffset, **conf)
        # load chromatic leakage
        if add_cl_vort is True:
            if verbose is True:
                print('   add chromatic leakage at the vortex plane')
            wf_cl = wf._wfarr*np.sqrt(vc_chrom_leak)
        else:
            wf_cl = 0
        # apply vortex
        scale_psf = wf._wfarr[0,0]/conf['psf_num'][0,0]
        wf_corr = (conf['psf_num']*conf['vvc'] - conf['perf_num'])*scale_psf
        wf._wfarr = wf._wfarr*conf['vvc'] - wf_corr + wf_cl
        # vortex phase mask physical size
        proper.prop_circular_aperture(wf, diam_vpm/2, 0, 0, NORM=False)
        # propagate to lyot stop (lt = light trap)
        lens(wf, zoffset_after=-vc_zoffset, 
            dist_start=lt_dist_start, diam_start=lt_diam_start,
            dist_end=lt_dist_end, diam_end=lt_diam_end, **conf)
    
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
        # lyot occulting mask physical size
        proper.prop_circular_aperture(wf, diam_lom/2, 0, 0, NORM=False)
        # propagate to lyot stop
        lens(wf, dist_start=lt_dist_start, diam_start=lt_diam_start,
            dist_end=lt_dist_end, diam_end=lt_diam_end, **conf)
    
    return wf