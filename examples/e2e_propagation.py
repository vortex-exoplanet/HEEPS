#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
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
import multiprocessing as mpro
from functools import partial
from sys import platform
import os

""" user inputs """
bands = ['Lp', 'Mp', 'N1', 'N2']
conf['CL_DIAM'] = 4
conf['ONAXIS'] = True
conf['STATIC_NCPA'] = False
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_subject'] = 'fenrir'
conf['send_message'] = 'HEEPS simulation finished'

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = None

# band specifications
band_specs = {'Lp': {'lam': 3.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'VC', 'RAVC']},
              'Mp': {'lam': 4.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'VC', 'RAVC']},
              'N1': {'lam': 8.7e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'VC']},
              'N2': {'lam': 11.5e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'VC']}}

# tip/tilt values
conf['tip_tilt'] = (0, 0)

# AO residual values
if True:
    AO_residuals_cube = fits.getdata(os.path.join(conf['INPUT_DIR'], conf['ATM_SCREEN_CUBE']))[:5]
#    AO_residuals_cube = fits.getdata('/mnt/disk4tb/METIS/SCAO/cube_compass_20181008_3600s_100ms.fits')[:9000]#[::3,:] # 300 ms sampling instead of 100 ms
else:
    AO_residuals_cube = np.array([None])
ncube = AO_residuals_cube.shape[0]
print('\nCube size = %s.'%ncube)

""" Create a function to propagate one single wavefront """

def propagate(wf_start, AO_residuals_cube, conf, ind):
    wf = copy.copy(wf_start)
    wavefront_abberations(wf, AO_residuals=AO_residuals_cube[ind], **conf)
    # METIS coronagraph modes
    if conf['MODE'] == 'ELT': # no Lyot stop (1, -0.3, 0)
        conf['LS_PARAMS'] = [1., -0.3, 0.]
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf, get_pupil='amp')
    elif conf['MODE'] in ('CL', 'CL4', 'CL5'): # classical lyot
        conf['LS_PARAMS'] = [0.8, 0.1, 1.1]
        if conf['ONAXIS'] == True:
            lyot(wf, conf)
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf)
    elif conf['MODE'] == 'VC':
        if conf['ONAXIS'] == True:
            vortex(wf, conf)
        #_, PUP = lyotstop(wf, conf, get_pupil='amp')
        lyotstop(wf, conf)
    elif conf['MODE'] == 'RAVC':
        RAVC = True
        apodization(wf, conf, RAVC=RAVC)
        if conf['ONAXIS'] == True:
            vortex(wf, conf)
        #_, PUP = lyotstop(wf, conf, RAVC=RAVC, get_pupil='amp')
        lyotstop(wf, conf, RAVC=RAVC)
    elif conf['MODE'] == 'APP': 
        APP = True
        if conf['ONAXIS'] == True:
            #_, PUP = lyotstop(wf, conf, APP=APP, get_pupil='phase')
            lyotstop(wf, conf, APP=APP)
    # get science image
    psf = detector(wf, conf)
    
    return psf

""" Start looping on the different bands """

for band in bands:
    
    # compute beam ratio, pupil size, and create the entrance pupil
    conf['WAVELENGTH'] = band_specs[band]['lam']
    conf['PIXEL_SCALE'] = band_specs[band]['pscale']
    wf_start, PUP = pupil(conf, get_pupil='amp')
    
    for i, mode in enumerate(band_specs[band]['modes']):
        
        # starting time
        t0 = time.time()
        
        # add mode to conf
        conf['MODE'] = mode
        
        # propagate the cube, using multiple cores if possible 
        if conf['cpucount'] != 1 and platform in ('linux', 'linux2', 'darwin'):
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
            psfs = np.zeros((ncube, conf['N_D'], conf['N_D']))
            for i in range(ncube):
                psfs[i,:,:] = propagate(wf_start, AO_residuals_cube, conf, i)
        
        # if only one frame, make dim = 2
        if ncube == 1:
            psfs = psfs[0]
        
        """ Write to .fits """
        
        # file name
        conf['PREFIX'] = ''
        filename = '%s%s_%s'%(conf['PREFIX'], conf['MODE'], band)
        # save PSF (or PSFs cube) and/or Lyot-Stop
        if True:
            fits.writeto(os.path.join(conf['OUT_DIR'], 'PSF_' + filename + '_4') \
                    + '.fits', np.float32(psfs), overwrite=True)
        if False:
            fits.writeto(os.path.join(conf['OUT_DIR'], 'LS_' + filename) \
                    + '.fits', PUP, overwrite=True)
        
        """ Save figures to .png """
        
        if False:
            plt.figure()
            #plt.imshow(psf**0.05, origin='lower')
            plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
            plt.colorbar()
            plt.show(block=False)
            plt.savefig(os.path.join(conf['OUT_DIR'], 'PSF_' + filename) \
                    + '.png', dpi=300, transparent=True)
        if False:
            plt.figure()
            plt.imshow(PUP[50:-50,50:-50], origin='lower', cmap='gray', vmin=0, vmax=1)
            #plt.imshow(PUP[50:-50,50:-50], origin='lower')
            plt.colorbar()
            plt.show(block=False)
            plt.savefig(os.path.join(conf['OUT_DIR'], 'LS_' + filename) \
                    + '.png', dpi=300, transparent=True)
        
        # print elapsed time
        print('      Elapsed %.3f seconds.'%(time.time() - t0))

# Send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: Simulation finished OK.\n", time.localtime()))
os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
        conf['send_subject'], conf['send_to']))
