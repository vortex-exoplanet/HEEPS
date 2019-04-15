#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import heeps
from heeps.config import conf
from heeps import metis_hci
import os.path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 


def sinusoid_ncpa_scaling(scale, amp, period, steps):
    '''Generates sinusoidal scaling values over a period of 1hr (60 min)
    with amplitude as amp and given period'''
    
    t_time = 60
    p = t_time/period
    start, end = -scale*amp, scale*amp
    time = np.linspace(0, t_time, steps)
    diff = end - start
    y = np.sin(2*np.pi*time/(p))*diff/2
    y += diff/2 + start
    return y

"""
Default simulation configuration defined in "read_config.py" can be 
overridden here by updating the dictionary
"""
conf['lam'] = 3.8e-6
conf['band'] = 'L'
conf['mode'] = 'CVC'
conf['prefix'] = 'test_'
conf['VC_charge'] = 2   # vortex charge is modified here
conf['onaxis'] = True   # True = on-axis, False = off-axis
conf['tip_tilt'] = [0, 0]
conf['petal_piston'] = [0,0,0,0,0,0]
conf['static_ncpa'] = True
conf['polish_error'] = False
conf['cpucount'] = 2    # specify CPU count

# loading cube of atmosphere phase screens
atm_screen_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))[0:5]

# Quasi-static linear NCPA creation over the length of atm screens
scale_val = 0.0089
ncpa_scaling_linear = np.linspace(-scale_val, scale_val, atm_screen_cube.shape[0])
# Quasi-static sine-wave NCPA creation over the length of atm screens
ncpa_sine_scaling = sinusoid_ncpa_scaling(scale_val, 3, 10, atm_screen_cube.shape[0]) 

# start end-to-end propagation of 
psf_cube = metis_hci(atm_screen_cube, ncpa_scaling_linear, **conf)
psf = psf_cube[0]

# figures
psf_filename = '%sPSF_%s_%s'%(conf['prefix'], conf['band'], conf['mode'])
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['output_dir'], psf_filename + '.png'), \
        dpi=300, transparent=True)

# save cube as fits file
fits.writeto(os.path.join(conf['output_dir'], psf_filename + '.fits'), \
        psf_cube, overwrite=True)
