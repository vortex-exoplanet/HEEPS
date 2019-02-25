#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

# Required libraries
import heeps
from heeps.config import conf
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
conf['tip_tilt'] = [0, 0]
conf['petal_piston'] = [0,0,0,0,0,0]
conf['static_ncpa'] = False
conf['polish_error'] = False

# loading one single atmosphere phase screen
atm_screen = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))[0]

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
heeps.aberrations.wavefront_aberrations(wfo, atm_screen=atm_screen, **conf)

#   3. Coronagraph selection -- e.g. RAVC, CVC, APP, CLC --
#    "Three coronagraphic planes: apodizer, focal plane mask & Lyot-stop"
heeps.coronagraphs.metis_modes(wfo, conf)

#   4. Detector plane
psf = heeps.detector.detector(wfo, conf)


""" Figures """
psf_filename = '%sPSF_%s_%s'%(conf['prefix'], conf['band'], conf['mode'])
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['output_dir'], psf_filename + '.png'), \
        dpi=300, transparent=True)
