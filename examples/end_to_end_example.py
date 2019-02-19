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


# HCI mode inputs
conf['lam'] = 3.8e-6
conf['band'] = 'L'
conf['mode'] = 'RAVC'
conf['VC_charge'] = 2 # vortex charge is modified here
conf['onaxis'] = True # True = on-axis, False = off-axis
conf['tip_tilt'] = [0., 0.]
conf['static_ncpa'] = False
conf['prefix'] = 'test_'

# loading cube of atmosphere phase screens
download_from_gdrive(conf['atm_screen_gdriveID'], conf['input_dir'], conf['atm_screen_cube'])
atm_screen = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_cube']))[0]

# =============================================================================
# End to end simulation example, showing propagation through each plane; 
#   1. Pupil creation 
#   2. Addition of wavefront aberrations
#   3. Coronagraph system
#   4. Detector plane
# =============================================================================

#   1. ELT Pupil Plane
wfo = heeps.pupil.pupil(conf)

#   2. Wavefront aberrations
heeps.aberrations.wavefront_aberrations(wfo, AO_residuals=atm_screen, **conf)

#   3. Coronagraph selection -- e.g. RAVC, CVC, APP, CLC --
#    "Three coronagraphic planes: apodizer, focal plane mask & Lyot-stop"
heeps.coronagraphs.metis_modes(wfo, conf)

#   4. Detector plane
psf = heeps.detector.detector(wfo, conf)


""" Figures """
filename_PSF = '%sPSF_%s_%s'%(conf['prefix'], conf['band'], conf['mode'])
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['output_dir'], filename_PSF + '.png'), \
        dpi=300, transparent=True)
