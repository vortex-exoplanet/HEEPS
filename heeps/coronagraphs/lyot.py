import heeps.util.img_processing as impro
import numpy as np
import proper
import os.path
from astropy.io import fits


def lyot(wfo, conf):
    
    f_lens = conf['focal']
    beam_ratio = conf['beam_ratio']
    CLC_diam = conf['CLC_diam'] # classical lyot diam in lam/D (default to 4)
    gridsize = conf['gridsize']
    tmp_dir = conf['temp_dir']
    
    proper.prop_propagate(wfo, f_lens, 'inizio') # propagate wavefront
    proper.prop_lens(wfo, f_lens, 'focusing lens LC') # apply lens
    proper.prop_propagate(wfo, f_lens, 'LC') # propagate wavefront
    
    # create or load the classical Lyot mask
    calib = str(CLC_diam)+str('_')+str(int(beam_ratio*100))+str('_')+str(gridsize)
    my_file = os.path.join(tmp_dir, 'clc_'+calib+'.fits')
    if not os.path.isfile(my_file):
        # calculate exact size of Lyot mask diameter, in pixels
        Dmask = CLC_diam/beam_ratio
        # oversample the Lyot mask (round up)
        samp = 100
        ndisk = int(samp*np.ceil(Dmask))
        ndisk = ndisk + 1 if not ndisk % 2 else ndisk # must be odd
        # find center
        cdisk = int((ndisk - 1)/2)
        # calculate the distances to center
        xy = range(-cdisk, cdisk + 1)
        x,y = np.meshgrid(xy, xy)
        dist = np.sqrt(x**2 + y**2)
        # create the Lyot mask
        mask = np.zeros((ndisk, ndisk))
        mask[np.where(dist > samp*Dmask/2)] = 1
        # resize to Lyot mask real size, and pad with ones
        mask = impro.resize_img(mask, int(ndisk/samp))
        mask = impro.pad_img(mask, gridsize, 1)
        # write mask
        fits.writeto(my_file, mask)
    else:
        mask = fits.getdata(my_file)
    
    # apply lyot mask
    mask = proper.prop_shift_center(mask)
    wfo._wfarr.real *= mask
    
    proper.prop_propagate(wfo, f_lens, "propagate to pupil reimaging lens")  
    proper.prop_lens(wfo, f_lens, "apply pupil reimaging lens")
    proper.prop_propagate(wfo, f_lens, "lyot stop")
    
    return wfo
