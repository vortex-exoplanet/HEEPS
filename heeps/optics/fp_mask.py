from .lens import lens
from .vortex_init import vortex_init
import os.path
import numpy as np
from heeps.util.img_processing import resize_img, pad_img
from astropy.io import fits
import proper

def fp_mask(wf, mode='RAVC', focal=660, verbose=False, **conf):

    # case 1: vortex coronagraphs
    if mode in ['CVC', 'RAVC']:
        if verbose is True:
            print('   apply vortex phase mask')                        
        # update conf
        conf.update(focal=focal)
        # load vortex calibration files: psf_num, vvc, perf_num
        conf = vortex_init(verbose=verbose, **conf)
        # propagate to vortex
        lens(wf, focal)
        # apply vortex
        scale_psf = wf._wfarr[0,0]/conf['psf_num'][0,0]
        wf_corr = (conf['psf_num']*conf['vvc'] - conf['perf_num'])*scale_psf
        wf._wfarr = wf._wfarr*conf['vvc'] - wf_corr
        # propagate to lyot stop
        lens(wf, focal)
    
    # case 2: classical Lyot
    elif mode in ['CLC']:
        if verbose is True:
            print('   apply classical lyot mask')
        # propagate to lyot mask
        lens(wf, focal)
        # create or load the classical Lyot mask
        calib = str(conf['clc_diam'])+str('_')+str(int(conf['beam_ratio']*100))+str('_')+str(conf['ngrid'])
        my_file = os.path.join(conf['dir_temp'], 'clc_'+calib+'.fits')
        if not os.path.isfile(my_file):
            # calculate exact size of Lyot mask diameter, in pixels
            Dmask = conf['clc_diam']/conf['beam_ratio']
            # oversample the Lyot mask (round up)
            samp = 100
            ndisk = int(samp*np.ceil(Dmask))
            ndisk = ndisk if ndisk % 2 else ndisk + 1# must be odd
            # find center
            cdisk = int((ndisk - 1)/2)
            # calculate the distances to center
            xy = range(-cdisk, cdisk + 1)
            x,y = np.meshgrid(xy, xy)
            r = abs(x + 1j*y)
            # create the Lyot mask
            mask = np.zeros((ndisk, ndisk))
            mask[np.where(r > samp*Dmask/2)] = 1
            # resize to Lyot mask real size, and pad with ones
            mask = resize_img(mask, int(ndisk/samp))
            mask = pad_img(mask, conf['ngrid'], 1)
            # write mask
            fits.writeto(my_file, mask, overwrite=True)
        else:
            mask = fits.getdata(my_file)
        
        # apply lyot mask
        mask = proper.prop_shift_center(mask)
        wf._wfarr.real *= mask
        
        # propagate to lyot stop
        lens(wf, focal)

    return wf