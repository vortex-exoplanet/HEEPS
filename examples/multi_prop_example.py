#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
import heeps
from heeps import metis_hci
from heeps.config import conf, download_from_gdrive
import os.path


"""
Default simulation configuration defined in "read_config.py" can be 
overridden here by updating the dictionary
"""
# HCI mode inputs
conf['lam'] = 3.8e-6
conf['band'] = 'L'
conf['mode'] = 'CVC'
conf['VC_charge'] = 2 # vortex charge is modified here
conf['onaxis'] = True # True = on-axis, False = off-axis
conf['tip_tilt'] = [0., 0.]
conf['static_ncpa'] = False
conf['prefix'] = 'test_'

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = 4

# your email here, to receive notification when finished
conf['send_to'] = None

# loading cube of atmosphere phase screens
download_from_gdrive(conf['atm_screen_gdriveID'], conf['input_dir'], conf['atm_screen_cube'])
atm_screen = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_cube']))[:5]

# start end-to-end propagation
psf = metis_hci(atm_screen, **conf)


""" Figures """
filename_PSF = '%sPSF_%s_%s'%(conf['prefix'], conf['band'], conf['mode'])
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['output_dir'], filename_PSF + '.png'), \
        dpi=300, transparent=True)
