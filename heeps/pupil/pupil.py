import heeps.util.img_processing as impro
import proper
import os.path
import numpy as np
import astropy.units as u
from astropy.io import fits

def pupil(conf, pupil_file=None, margin=0):
    
    # load useful parameters
    lam = conf['lam']
    diam = conf['diam']
    gridsize = conf['gridsize']
    beam_ratio = (conf['pscale']*u.mas/(lam/diam)).to('rad').value # lam/D per pix
    get_amp = conf['get_amp']
    get_phase = conf['get_phase']
    
    # begin PROPER
    wf = proper.prop_begin(diam, lam, gridsize, beam_ratio)
    
    # store the beam ratio and the pupil size in conf
    # pupil size must be odd (PROPER sets the center up-right next to the grid center)
    conf['beam_ratio'] = beam_ratio
    npupil = np.ceil(gridsize*beam_ratio)
    npupil = int(npupil + 1) if npupil % 2 == 0 else int(npupil)
    conf['npupil'] = npupil
    
    # select a pupil file, amongst various missing segments configurations
    pupil_file = {
    0: conf['pupil_file'],
    1: 'ELT_2048_37m_11m_5mas_nospiders_1missing_cut.fits',
    2: 'ELT_2048_37m_11m_5mas_nospiders_2missing_cut.fits',
    4: 'ELT_2048_37m_11m_5mas_nospiders_4missing_cut.fits',
    7: 'ELT_2048_37m_11m_5mas_nospiders_7missing_cut_flower.fits'}
    
    # create a pupil with circular obscuration (default), or load pupil from file
    if conf['pupil_file'] is None and conf['N_mis_segments'] == 0:
        proper.prop_circular_aperture(wf, diam/2)
        proper.prop_circular_obscuration(wf, conf['R_obstr'], NORM=True)
    else:
        # load pupil, resize, and pad with zeros to match PROPER gridsize
        pupil = fits.getdata(os.path.join(conf['input_dir'], pupil_file[conf['N_mis_segments']]))
        pupil = impro.resize_img(pupil, npupil)
        pupil = impro.pad_img(pupil, gridsize)
        # multiply the loaded pupil
        proper.prop_multiply(wf, pupil)
        
    # add spiders
    if (conf['spiders_width'] != 0):
        for angle in conf['spiders_angle']:
            proper.prop_rectangular_obscuration(wf, conf['spiders_width'], 2*diam, ROTATION=angle)
    
    # define the entrance wavefront
    proper.prop_define_entrance(wf)
    wf.wfarr /= np.amax(wf._wfarr) # max(amplitude)=1
    
    # get the pupil amplitude and phase for output
    pup_amp = impro.crop_img(proper.prop_get_amplitude(wf), npupil, margin)\
            if get_amp is True else None
    pup_phase = impro.crop_img(proper.prop_get_phase(wf), npupil, margin)\
            if get_phase is True else None
    
    return wf, pup_amp, pup_phase
