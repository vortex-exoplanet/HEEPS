#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
from read_config import conf
from heeps import metis_hci, download_from_gdrive

"""
Default simulation cofiguration defined in "read_config.py" can be 
overridden here by updating the dictionary
"""
conf
conf['WAVELENGTH'] = 3.80*10**-6 
conf['CHARGE'] = 2 # charge is modified here


#atm_screen = 0 # no atmospheric phase screen

# 2-D pahse screen
#atm_screen = fits.getdata(conf['INPUT_DIR'] + conf['ATM_SCREEN_2D']) 

# getting multi-cube phase screen
download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
atm_screen = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0:10] 

TILT = np.array(conf['TILT_2D']) 
#TILT = np.random.randn(conf['TILT_CUBE'],2)


psf = metis_hci(conf['MODE'], conf, atm_screen, TILT)            


plt.figure()
plt.imshow(psf**0.5)
plt.colorbar()

plt.show()

