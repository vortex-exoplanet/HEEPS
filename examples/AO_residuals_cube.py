#!/usr/bin/env python3

# =============================================================================
#       Example file for creating a coronagraphic/non-coronagraphic METIS PSF
# =============================================================================

""" Required libraries """
import matplotlib.pyplot as plt # for plotting simulated PSFs
import numpy as np
from astropy.io import fits 
from heeps.config import conf, download_from_gdrive
from heeps.pupil import pupil
from heeps.abberations import wavefront_abberations
from heeps.coronagraphs import apodization, vortex, lyotstop, lyot
from heeps.detector import detector
import os.path
import copy
import time


""" Default simulation configuration defined in "read_config.py" can be 
overridden here by updating the dictionary """
conf['WAVELENGTH'] = 3.80*10**-6 
conf['STATIC_NCPA'] = False
conf['MODE'] = 'VC'
conf['CL_DIAM'] = 4 # classical lyot diam in lam/D

""" Pupil """
conf['PUPIL_FILE'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
wfo, PUP = pupil(conf, get_pupil='amp')

""" Loading tip-tilt values """
conf['tip_tilt'] = (0, 0)
# conf['tip_tilt'] = np.random.randn(conf['TILT_CUBE'], 2)

""" Loading AO residual values, getting multi-cube phase screen from Google Drive """
if True:
    download_from_gdrive(conf['GDRIVE_ID'], conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE'])
    AO_residuals_cube = fits.getdata(conf['INPUT_DIR'] + conf['ATM_SCREEN_CUBE'])[1:15]
else:
    AO_residuals_cube = np.array([None])
ncube = AO_residuals_cube.shape[0]
print('ncube = %s'%ncube)

psfs = None
t0 = time.time()
for i in range(ncube):
    print(i)
    wf = copy.copy(wfo)
    wavefront_abberations(wf, AO_residuals=AO_residuals_cube[i], **conf)
    
    """ Coronagraph """
    if conf['MODE'] == 'ELT': # no Lyot stop (1, -0.3, 0)
        conf['LS_PARAMS'] = [1., -0.3, 0.]
        lyotstop(wf, conf)
    elif conf['MODE'] == 'OFFAXIS': # only Lyot stop
        lyotstop(wf, conf)
    elif conf['MODE'] == 'OFFAXIS_RA': # only Lyot stop and ring apodizer
        RAVC=True
        apodization(wf, conf, RAVC=RAVC)
        lyotstop(wf, conf, RAVC=RAVC)
    elif conf['MODE'] == 'CL': # classical lyot
        conf['LS_PARAMS'] = [0.8, 0.1, 1.1]
        lyot(wf, conf)
        lyotstop(wf, conf, RAVC=False)
    elif conf['MODE'] == 'VC':
        vortex(wf, conf)
        lyotstop(wf, conf, RAVC=False)
    elif conf['MODE'] == 'RAVC':
        RAVC=True
        apodization(wf, conf, RAVC=RAVC)
        vortex(wf, conf)
        lyotstop(wf, conf, RAVC=RAVC)
    elif conf['MODE'] == 'APP': 
        _, PUP = lyotstop(wf, conf, APP=True, get_pupil='phase')
    
    """ Science image """
    psf = detector(wf, conf)
    
    # stack psfs
    if psfs is None:
        psfs = psf[None, ...]
    else:
        psfs = np.concatenate([psfs, psf[None, ...]])
print(time.time() - t0)

""" write cube to fits """
if True:
    filename_PSF_cube = 'NEW_compass_10min_100ms_' + conf['MODE']
    fits.writeto(os.path.join(conf['OUT_DIR'], filename_PSF_cube) + '.fits', psfs, overwrite=True)



""" write to fits """
if False:
    conf['PREFIX'] = ''
    filename_PSF = conf['PREFIX'] + conf['MODE'] + '_PSF'
    filename_LS = conf['PREFIX'] + conf['MODE'] + '_LS'
    fits.writeto(os.path.join(conf['OUT_DIR'], filename_PSF) + '.fits', psf, overwrite=True)
    fits.writeto(os.path.join(conf['OUT_DIR'], filename_LS) + '.fits', PUP, overwrite=True)

""" Figures """
if False:
    plt.figure(1)
    #plt.imshow(psf**0.05, origin='lower')
    plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
    plt.colorbar()
    plt.show(block=False)
    plt.savefig(os.path.join(conf['OUT_DIR'], filename_PSF) + '.png', dpi=300, transparent=True)

    plt.figure(2)
    plt.imshow(PUP[50:-50,50:-50], origin='lower', cmap='gray', vmin=0, vmax=1)
    #plt.imshow(PUP[50:-50,50:-50], origin='lower')
    plt.colorbar()
    plt.show(block=False)
    plt.savefig(os.path.join(conf['OUT_DIR'], filename_LS) + '.png', dpi=300, transparent=True)
