#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
from heeps import metis_hci
from heeps.config import conf, download_from_gdrive


"""
Default simulation cofiguration defined in "read_config.py" can be 
overridden here by updating the dictionary
"""
conf['WAVELENGTH'] = 3.80*10**-6 
conf['MODE'] = 'VC'

# getting multi-cube phase screen

get_cube = True
if get_cube is True:
    download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
    conf['ATM_SCREEN'] = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0:5] 

TILT = np.array(conf['TILT_2D']) 

psf = metis_hci(conf, conf['ATM_SCREEN']  , TILT)


plt.figure()
#plt.imshow(psf**0.25, origin='lower')
plt.imshow(np.log10(psf), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show()


