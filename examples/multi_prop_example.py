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
conf['tip_tilt'] = [0., 0.]
conf['static_ncpa'] = False
conf['cpucount'] = 8    # specify CPU count

# loading cube of atmosphere phase screens
atm_screen_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))

# start end-to-end propagation of 
psf_cube = metis_hci(atm_screen_cube, **conf)
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
