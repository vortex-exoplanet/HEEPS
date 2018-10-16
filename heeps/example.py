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
from heeps.coronagraphs import apodization, vortex, lyotstop, lyot
from heeps.detector import detector
import os.path

""" Default simulation configuration defined in "read_config.py" can be 
overridden here by updating the dictionary """
conf['WAVELENGTH'] = 3.80*10**-6 
conf['STATIC_NCPA'] = False
conf['TILT_2D'] = np.zeros(2) #np.random.randn(conf['TILT_CUBE'],2)
conf['ATM_SCREEN'] = 0
conf['PHASE_APODIZER_FILE'] = 0
conf['MODE'] = 'ELT'
conf['DIAM_CL'] = 4 # classical lyot diam in lam/D
pre_sim(conf) 

""" Pupil """
conf['PUPIL_FILE'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
wfo = pupil(conf) 

""" Wavefront errors, getting multi-cube phase screen from Google Drive"""
get_cube = True
if get_cube is True:
    download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']) 
    conf['ATM_SCREEN'] = fits.getdata(conf['INPUT_DIR']+ conf['ATM_SCREEN_CUBE'])[0] 
wfo = wavefront_abberations(wfo, conf, conf['ATM_SCREEN'], conf['TILT_2D'])

""" Coronagraph """
if conf['MODE'] == 'ELT': # no Lyot stop (1, -0.3, 0)
    conf['LS_PARA'] = [1., -0.3, 0.]
    _, LS_pupil = lyotstop(wfo, conf, RAVC=False)
elif conf['MODE'] == 'OFFAXIS': # only lyot stop
    lyotstop(wfo, conf, RAVC=False)
elif conf['MODE'] == 'CL': # classical lyot
    lyot(wfo, conf)
    _, LS_pupil = lyotstop(wfo, conf, RAVC=False)
elif conf['MODE'] == 'VC':
    vortex(wfo, conf)
    lyotstop(wfo, conf, RAVC=False)
elif conf['MODE'] == 'RAVC':
    apodization(wfo, conf, RAVC=True)
    vortex(wfo, conf)
    lyotstop(wfo, conf, RAVC=True)
elif conf['MODE'] == 'MASK': # only ring apodizer and lyot stop
    apodization(wfo, conf, RAVC=True)
    lyotstop(wfo, conf, RAVC=True)
elif conf['MODE'] == 'APP': 
    conf['PHASE_APODIZER_FILE'] = 'app_phase_cut.fits'
    lyotstop(wfo, conf, RAVC=False)

""" Science image """
filename_PSF = conf['PREFIX'] + '_PSF_' + conf['MODE']
filename_LS = conf['PREFIX'] + '_LS_' + conf['MODE']
psf = detector(wfo, conf)
fits.writeto(os.path.join(conf['OUT_DIR'], filename_PSF) + '.fits', psf, overwrite=True)

""" Figures """
plt.figure(1)
#plt.imshow(psf**0.05, origin='lower')
plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['OUT_DIR'], filename_PSF) + '.png', dpi=300, transparent=True)

LS_pupil = LS_pupil[50:-50,50:-50]
plt.figure(2)
plt.imshow(LS_pupil, origin='lower', cmap='gray', vmin=0, vmax=1)
plt.colorbar()
plt.show(block=False)
plt.savefig(os.path.join(conf['OUT_DIR'], filename_LS) + '.png', dpi=300, transparent=True)
