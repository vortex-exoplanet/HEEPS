#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

""" Required libraries """
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
from heeps.config import conf, pre_sim, download_from_gdrive
from heeps.pupil import pupil
from heeps.abberations import wavefront_abberations
from heeps.coronagraphs import apodization, vortex, lyotstop
from heeps.detector import detector

""" Default simulation configuration defined in "read_config.py" can be 
overridden here by updating the dictionary """
conf['WAVELENGTH'] = 3.80*10**-6 
conf['STATIC_NCPA'] = False
conf['TILT_2D'] = np.zeros(2) #np.random.randn(conf['TILT_CUBE'],2)
conf['ATM_SCREEN'] = 0
conf['RAVC'] = False
conf['PHASE_APODIZER_FILE'] = 0
pre_sim(conf) 

""" Pupil """
wfo = pupil(conf) 

""" Wavefront errors, getting multi-cube phase screen from Google Drive"""
get_cube = True
if get_cube is True:
	download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
	conf['ATM_SCREEN'] = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0] 
wfo = wavefront_abberations(wfo, conf, conf['ATM_SCREEN'], conf['TILT_2D'])

""" Coronagraph """
RAVC = False
wfo = apodization(wfo, conf, conf['RAVC'])
wfo = vortex(wfo, conf)
wfo = lyotstop(wfo, conf, conf['RAVC'])

""" Science image """
psf = detector(wfo, conf)
fits.writeto(conf['OUT_DIR'] + conf['PREFIX'] + '_PSF_'+ conf['MODE'] +'.fits', psf, overwrite=True)        

""" Figure """
plt.figure()
plt.imshow(psf**0.05)
plt.colorbar()
plt.show(block=False)

