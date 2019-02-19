import numpy as np
import proper
from skimage.transform import resize
from astropy.io import fits
import os.path

def lyotstop(wf, conf, RAVC=None, APP=None, get_pupil='no', dnpup=50):
    """Add a Lyot stop, or an APP."""
    
    # load parameters
    npupil = conf['npupil']
    pad = int((conf['gridsize'] - npupil)/2)
    
    # get LS misalignments
    LS_misalignment = (np.array(conf['LS_misalign'])*npupil).astype(int)
    dx_amp, dy_amp, dz_amp = LS_misalignment[0:3]
    dx_phase, dy_phase, dz_phase = LS_misalignment[3:6]
    
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
        APP_amp = resize(APP_amp, (npupil, npupil), preserve_range=True, mode='reflect')
        APP_phase = resize(APP_phase, (npupil, npupil), preserve_range=True, mode='reflect')
        # pad with zeros to match PROPER gridsize
        APP_amp = np.pad(APP_amp, [(pad+1+dx_amp, pad-dx_amp), \
                (pad+1+dy_amp, pad-dy_amp)], mode='constant')
        APP_phase = np.pad(APP_phase, [(pad+1+dx_phase, pad-dx_phase), \
                (pad+1+dy_phase, pad-dy_phase)], mode='constant')
        # multiply the loaded APP
        proper.prop_multiply(wf, APP_amp*np.exp(1j*APP_phase))
    
    # get the pupil amplitude or phase for output
    if get_pupil.lower() in 'amplitude':
        return wf, proper.prop_get_amplitude(wf)[pad+1-dnpup:-pad+dnpup, pad+1-dnpup:-pad+dnpup]
    elif get_pupil.lower() in 'phase':
        return wf, proper.prop_get_phase(wf)[pad+1-dnpup:-pad+dnpup, pad+1-dnpup:-pad+dnpup]
    else:
        return wf
