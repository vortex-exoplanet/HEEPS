#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 

from heeps.config import download_from_gdrive
from heeps import metis_hci
from heeps.config import conf, pre_sim

pre_sim(conf) 

"""
Default simulation cofiguration defined in "read_config.py" can be 
overridden here by updating the dictionary
"""
conf['WAVELENGTH'] = 3.80*10**-6 
conf['CHARGE'] = 2 # charge is modified here
conf['MODE'] = 'VC'
conf['STATIC_NCPA'] = False

# 2-D pahse screen
atm_screen = 0 

# getting multi-cube phase screen
download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
atm_screen = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0] 

TILT = np.array(conf['TILT_2D']) 
#TILT = np.random.randn(conf['TILT_CUBE'],2)

psf = metis_hci(conf, atm_screen, TILT)


plt.figure()
plt.imshow(psf**0.05)
plt.colorbar()

plt.show()

