#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits 
from heeps.config import conf, download_from_gdrive
from heeps.pupil import pupil
from heeps.aberrations import wavefront_aberrations
from heeps.coronagraphs import apodization, vortex, lyotstop, lyot
from heeps.detector import detector
import os.path
from copy import deepcopy
import time
import multiprocessing as mpro
from functools import partial
from sys import platform
import os

""" user inputs """
conf['bands'] = ['L', 'M', 'N1', 'N2']
conf['CLC_diam'] = 4
conf['onaxis'] = True
conf['static_ncpa'] = False
conf['send_to'] = None #'cdelacroix@uliege.be'
conf['send_message'] = 'HEEPS simulation finished.'

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = None

# band specifications
band_specs = {'L': {'lam': 3.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'CVC', 'RAVC']},
              'M': {'lam': 4.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'CVC', 'RAVC']},
              'N1': {'lam': 8.7e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'CVC']},
              'N2': {'lam': 11.5e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'CVC']}}

# tip/tilt values
conf['tip_tilt'] = (0, 0)

# AO residual values
if True:
    AO_residuals_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_cube']))[:5]
#    AO_residuals_cube = fits.getdata('/mnt/disk4tb/METIS/SCAO/cube_compass_20181008_3600s_100ms.fits')[::3,:]
else:
    AO_residuals_cube = np.array([None])
ncube = AO_residuals_cube.shape[0]
print('\nCube size = %s.'%ncube)

""" Create a function to propagate one single wavefront """

def propagate(wf_start, AO_residuals_cube, conf, ind):
    wf = deepcopy(wf_start)
    wavefront_aberrations(wf, AO_residuals=AO_residuals_cube[ind], **conf)
    # METIS coronagraph modes
    if conf['mode'] == 'ELT': # no Lyot stop (1, -0.3, 0)
        conf['LS_params'] = [1., -0.3, 0.]
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf, get_pupil='amp')
    elif conf['mode'] == 'CLC': # classical lyot
        conf['LS_params'] = [0.8, 0.1, 1.1]
        if conf['onaxis'] == True:
            lyot(wf, conf)
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf)
    elif conf['mode'] == 'CVC':
        if conf['onaxis'] == True:
            vortex(wf, conf)
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf)
    elif conf['mode'] == 'RAVC':
        RAVC = True
        apodization(wf, conf, RAVC=RAVC)
        if conf['onaxis'] == True:
            vortex(wf, conf)
        #_, PUP = lyotstop(wf, conf, RAVC=RAVC, get_pupil='amp')
        lyotstop(wf, conf, RAVC=RAVC)
    elif conf['mode'] == 'APP': 
        APP = True
        if conf['onaxis'] == True:
            #_, PUP = lyotstop(wf, conf, APP=APP, get_pupil='phase')
            lyotstop(wf, conf, APP=APP)
    # get science image
    psf = detector(wf, conf)
    
    return psf

""" Start looping on the different bands """

for band in conf['bands']:
    conf['band'] = band
    conf['lam'] = band_specs[band]['lam']
    conf['pscale'] = band_specs[band]['pscale']
    
    # compute beam ratio, pupil size, and create the entrance pupil
    wf_start, PUP = pupil(conf, get_pupil='amp')
    
    for i, mode in enumerate(band_specs[band]['modes']):
        conf['mode'] = mode
        
        # starting time
        t0 = time.time()
        
        # propagate the cube, using multiple cores if possible 
        if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
            if conf['cpucount'] == None:
                conf['cpucount'] = mpro.cpu_count() - 1
            print('%s: %s band, %s mode, using %s cores.'\
                    %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                    band, mode, conf['cpucount']))
            p = mpro.Pool(conf['cpucount'])
            func = partial(propagate, wf_start, AO_residuals_cube, conf)
            psfs = np.array(p.map(func, range(ncube)))
        else:
            print('%s: %s band, %s mode, using %s core.'\
                    %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                    band, mode, 1))
            psfs = np.zeros((ncube, conf['ndet'], conf['ndet']))
            for i in range(ncube):
                psfs[i,:,:] = propagate(wf_start, AO_residuals_cube, conf, i)
        
        # if only one frame, make dim = 2
        if ncube == 1:
            psfs = psfs[0]
        
        """ Write to .fits """
        
        # file name
        conf['prefix'] = ''
        filename = '%s%s_%s'%(conf['prefix'], conf['mode'], band)
        # save PSF (or PSFs cube) and/or Lyot-Stop
        if True:
            fits.writeto(os.path.join(conf['output_dir'], 'PSF_' + filename) \
                    + '.fits', np.float32(psfs), overwrite=True)
        if False:
            fits.writeto(os.path.join(conf['output_dir'], 'LS_' + filename) \
                    + '.fits', PUP, overwrite=True)
        
        """ Save figures to .png """
        
        if False:
            plt.figure()
            #plt.imshow(psf**0.05, origin='lower')
            plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
            plt.colorbar()
            plt.show(block=False)
            plt.savefig(os.path.join(conf['output_dir'], 'PSF_' + filename) \
                    + '.png', dpi=300, transparent=True)
        if False:
            plt.figure()
            plt.imshow(PUP[50:-50,50:-50], origin='lower', cmap='gray', vmin=0, vmax=1)
            #plt.imshow(PUP[50:-50,50:-50], origin='lower')
            plt.colorbar()
            plt.show(block=False)
            plt.savefig(os.path.join(conf['output_dir'], 'LS_' + filename) \
                    + '.png', dpi=300, transparent=True)
        
        # print elapsed time
        print('      Elapsed %.3f seconds.'%(time.time() - t0))

# Send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: End-to-end simulation finished OK.\n", time.localtime()))
os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
        conf['send_subject'], conf['send_to']))
