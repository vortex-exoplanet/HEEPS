import heeps.util.img_processing as impro
import numpy as np
import proper
from astropy.io import fits
import os.path

def lyotstop(wf, conf, RAVC=None, APP=None, margin=50):
    """Add a Lyot stop, or an APP."""
    
    # load useful parameters
    npupil = conf['npupil']
    gridsize = conf['gridsize']
    get_amp = conf['get_amp']
    get_phase = conf['get_phase']
    
    # get LS misalignments
#    LS_misalign = np.int64(np.array(conf['LS_misalign'])*npupil)
    LS_misalign = conf['LS_misalign'] if conf['LS_misalign'] else [0,0,0,0,0,0]
    dx_amp, dy_amp, dz_amp = LS_misalign[0:3]
    dx_phase, dy_phase, dz_phase = LS_misalign[3:6]
    
    # case 1: Lyot stop (no APP)
    if APP is not True:
        
        # Lyot stop parameters: R_out, dR_in, spi_width
        # outer radius (absolute %), inner radius (relative %), spider width (m)
        (R_out, dR_in, spi_width) = conf['LS_params']
        
        # Lyot stop inner radius at least as large as obstruction radius
        R_in = conf['R_obstr']
        
        # case of a ring apodizer
        if RAVC is True:
            # define the apodizer transmission and apodizer radius [Mawet2013]
            # apodizer radius at least as large as obstruction radius
            T_ravc = 1 - (R_in**2 + R_in*np.sqrt(R_in**2 + 8))/4
            R_in /= np.sqrt(1 - T_ravc)
        
        # oversize Lyot stop inner radius
        R_in += dR_in
        
        # create Lyot stop
        proper.prop_circular_aperture(wf, R_out, dx_amp, dy_amp, NORM=True)
        if R_in > 0:
            proper.prop_circular_obscuration(wf, R_in, dx_amp, dy_amp, NORM=True)
        if spi_width > 0:
            for angle in conf['spiders_angle']:
                proper.prop_rectangular_obscuration(wf, spi_width, 2*conf['diam'], \
                        dx_amp, dy_amp, ROTATION=angle)
    
    # case 2: APP (no Lyot stop)
    else: 
        # get amplitude and phase files
        APP_amp_file = os.path.join(conf['input_dir'], conf['APP_amp_file'])
        APP_phase_file = os.path.join(conf['input_dir'], conf['APP_phase_file'])
        # get amplitude and phase data
        APP_amp = fits.getdata(APP_amp_file) if os.path.isfile(APP_amp_file) \
                else np.ones((npupil, npupil))
        APP_phase = fits.getdata(APP_phase_file) if os.path.isfile(APP_phase_file) \
                else np.zeros((npupil, npupil))
        # resize to npupil
        APP_amp = impro.resize_img(APP_amp, npupil)
        APP_phase = impro.resize_img(APP_phase, npupil)
        # pad with zeros to match PROPER gridsize
        APP_amp = impro.pad_img(APP_amp, gridsize, 1)
        APP_phase = impro.pad_img(APP_phase, gridsize, 0)
        # multiply the loaded APP
        proper.prop_multiply(wf, APP_amp*np.exp(1j*APP_phase))
    
    # get the LS amplitude and phase for output
    LS_amp = impro.crop_img(proper.prop_get_amplitude(wf), npupil, margin)\
            if get_amp is True else None
    LS_phase = impro.crop_img(proper.prop_get_phase(wf), npupil, margin)\
            if get_phase is True else None
    
    return wf, LS_amp, LS_phase
