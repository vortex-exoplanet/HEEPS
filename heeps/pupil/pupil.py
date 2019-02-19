import numpy as np
from skimage.transform import resize
import proper
from astropy.io import fits
import astropy.units as u
import os.path

def pupil(conf, pupil_file=None, get_pupil='no', dnpup=50):
    
    # load useful parameters
    lam = conf['lam']
    diam = conf['diam']
    gridsize = conf['gridsize']
    
    # compute the beam ratio
    beam_ratio = (conf['pscale']*u.mas/(lam/diam)).to('rad').value
    conf['beam_ratio'] = beam_ratio
    # compute the pupil size, must be odd (PROPER sets the center up-right next to the grid center)
    npupil = np.ceil(gridsize*beam_ratio)
    npupil = int(npupil + 1) if npupil % 2 == 0 else int(npupil)
    conf['npupil'] = npupil
    # calculate the size of the pupil to pad with zeros
    pad = int((gridsize - npupil)/2)
    
    # begin PROPER
    wfo = proper.prop_begin(diam, lam, gridsize, beam_ratio)
    
    # select pupil file, amongst various missing segments configurations
    pupil_file = {
    0: conf['pupil_file'],
    1: '/ELT_2048_37m_11m_5mas_nospiders_1missing_cut.fits',
    2: '/ELT_2048_37m_11m_5mas_nospiders_2missing_cut.fits',
    4: '/ELT_2048_37m_11m_5mas_nospiders_4missing_cut.fits',
    7: '/ELT_2048_37m_11m_5mas_nospiders_7missing_1_cut.fits'}
    
    # create a pupil with circular obscuration (default), or load pupil from file
    if pupil_file is None:
        proper.prop_circular_aperture(wfo, diam/2)
        proper.prop_circular_obscuration(wfo, conf['R_obstr'], NORM=True)
    else:
        # load pupil, resize, and pad with zeros to match PROPER gridsize
        pupil = fits.getdata(os.path.join(conf['input_dir'], pupil_file[conf['N_mis_segments']]))
        pupil = resize(pupil, (npupil, npupil), preserve_range=True, mode='reflect')
        pupil = np.pad(pupil, [(pad+1,pad),(pad+1,pad)], mode='constant')
        # multiply the loaded pupil
        proper.prop_multiply(wfo, pupil)
        
    # add spiders
    if (conf['spiders_width'] != 0):
        for angle in conf['spiders_angle']:
            proper.prop_rectangular_obscuration(wfo, conf['spiders_width'], 2*diam, ROTATION=angle)
    
    # define the entrance wavefront
    proper.prop_define_entrance(wfo)
    wfo.wfarr /= np.amax(wfo._wfarr) # max(amplitude)=1
    
    # get the pupil amplitude or phase for output
    if get_pupil.lower() in 'amplitude':
        return wfo, proper.prop_get_amplitude(wfo)[pad+1-dnpup:-pad+dnpup, pad+1-dnpup:-pad+dnpup]
    elif get_pupil.lower() in 'phase':
        return wfo, proper.prop_get_phase(wfo)[pad+1-dnpup:-pad+dnpup, pad+1-dnpup:-pad+dnpup]
    else:
        return wfo
