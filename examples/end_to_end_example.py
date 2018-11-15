#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 

import heeps
from heeps.config import conf, download_from_gdrive
import os.path


conf['WAVELENGTH'] = 3.80*10**-6 
conf['CHARGE'] = 2 # charge is modified here
conf['MODE'] = 'MASK'
conf['STATIC_NCPA'] = False


# loading multi-cube phase screen
get_cube = True
if get_cube is True:
    download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
    conf['ATM_SCREEN'] = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0] 

tilt = np.array(conf['TILT_2D']) 


# =============================================================================
# End to end simulation example, showing propagation through each plane; 
#   1. Pupil creation 
#   2. Addition of wavefront abberations
#   3. Coronagraph system
#   4. Detector plane
# =============================================================================


#   1. ELT Pupil Plane
wfo = heeps.pupil.pupil(conf)

#   2. Wavefront abberations
heeps.abberations.wavefront_abberations(wfo, AO_residuals=conf['ATM_SCREEN'], 
        tip_tilt=tilt, **conf)

#   3. Coronagraph selection -- Vortex Classical (VC) / RAVC / CL / APP --
#    "Three coronagraphic planes Apodizer, focal plane mask & Lyot-stop"
heeps.coronagraphs.metis_modes(wfo, conf)

#   4. Detector plane
psf = heeps.detector.detector(wfo,conf)	

""" Science image """
filename_PSF = conf['PREFIX'] + '_PSF_' + conf['MODE']


""" Figures """
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show()
plt.savefig(os.path.join(conf['OUT_DIR'], filename_PSF) + '.png', dpi=300, transparent=True)





