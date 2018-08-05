#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
from read_config import conf
from heeps import metis_hci, download_from_gdrive


conf

"""
Default simulation cofiguration defined in "simulation_config.md" can be 
overridden here by updating the dictionary
"""

conf['WAVELENGTH'] = 3.80*10**-6 
conf['CHARGE'] = 2 # charge is modified here

atm_screen = fits.getdata(conf['INPUT_DIR'] + conf['ATM_SCREEN_2D'])

gdrive_id = '1AUtELRfn_xjnbsMM_SJG6W0c26zyzqNH'    # [Required] google drive id linked to the fits file 
atm_screen_cube  = conf['ATM_SCREEN_CUBE']
download_from_gdrive(gdrive_id, conf['INPUT_DIR'], atm_screen_cube) 
atm_screen = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0:10] 

TILT = np.array(conf['TILT_2D'])
#TILT = np.random.randn(10,2)


mode = 'VC'
conf['MODE'] = mode

psf = metis_hci(mode , conf, atm_screen, TILT)            


plt.figure()
plt.imshow(psf**0.25)
plt.colorbar()

plt.show()

