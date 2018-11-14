import numpy as np
from skimage.transform import resize
import proper
from astropy.io import fits
import os.path

def pupil(conf, pupil_file=None):
    
    # load useful parameters
    lam = conf['WAVELENGTH']
    diam = conf['DIAM']
    gridsize = conf['GRIDSIZE']
    
    # compute the beam ratio
    beam_ratio = conf['PIXEL_SCALE']*4.85e-9/(lam/diam)
    conf['beam_ratio'] = beam_ratio
    # compute the pupil size, must be odd (PROPER sets the center up-right next to the grid center)
    npupil = np.ceil(gridsize*beam_ratio)
    npupil = int(npupil + 1) if npupil % 2 == 0 else int(npupil)
    conf['NPUPIL'] = npupil
    
    # select pupil file, amongst various missing segments configurations
    pupil_file = {
    0: conf['PUPIL_FILE'],
    1: '/ELT_2048_37m_11m_5mas_nospiders_1missing_cut.fits',
    2: '/ELT_2048_37m_11m_5mas_nospiders_2missing_cut.fits',
    4: '/ELT_2048_37m_11m_5mas_nospiders_4missing_cut.fits',
    7: '/ELT_2048_37m_11m_5mas_nospiders_7missing_1_cut.fits'}
    
    # begin PROPER
    wfo = proper.prop_begin(diam, lam, gridsize, beam_ratio)
    
    # create a pupil with circular obscuration (default), or load pupil from file
    if pupil_file is None:
        proper.prop_circular_aperture(wfo, diam/2)
        proper.prop_circular_obscuration(wfo, conf['R_OBSTR'], NORM=True)
    else:
        # load pupil, resize, and pad with zeros to match PROPER gridsize
        pupil = fits.getdata(os.path.join(conf['INPUT_DIR'], pupil_file[conf['MIS_SEGMENTS_NU']]))
        pupil = resize(pupil, (npupil, npupil), preserve_range=True, mode='reflect')
        r = int((proper.prop_get_gridsize(wfo) - npupil)/2)
        pupil = np.pad(pupil, [(r+1,r),(r+1,r)], mode='constant')
        # multiply the loaded pupil
        proper.prop_multiply(wfo, pupil)
        
    # add spiders
    if (conf['SPIDERS_WIDTH'] != 0):
        for angle in conf['SPIDERS_ANGLE']:
            proper.prop_rectangular_obscuration(wfo, conf['SPIDERS_WIDTH'], 2*diam, ROTATION=angle)
    
    # define the entrance wavefront
    proper.prop_define_entrance(wfo)
    wfo.wfarr /= np.amax(wfo._wfarr) # max(amplitude)=1
    
    return wfo
