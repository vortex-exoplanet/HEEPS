#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits 
from heeps.config import conf
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
from itertools import zip_longest
import astropy.units as u

""" user inputs """

conf['bands'] = ['L', 'M', 'N1', 'N2']
conf['CLC_diam'] = 4
conf['onaxis'] = True
conf['static_ncpa'] = False
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_message'] = 'End-to-end simulation finished OK.'
conf['get_amp'] = False
conf['get_phase'] = False
conf['APP_phase_file'] = 'rot90_METIS_vAPP_PDR_L.fits'

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = None

# band specifications
band_specs = {'L': {'lam': 3.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
              'M': {'lam': 4.8e-6,
                  'pscale': 5.21,
                   'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
              'N1': {'lam': 8.7e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'CVC', 'CLC']},
              'N2': {'lam': 11.5e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'CVC', 'CLC']}}

""" Create a function to propagate one single wavefront """

def propagate(wf_start, conf, full_output, atm_screen, ncpa_screen, 
            petal_piston, misalign, zernike):
    # keep a copy of the input pupil wavefront
    wf = deepcopy(wf_start)
    # apply wavefront aberrations
    wavefront_aberrations(wf, zernike=zernike, petal_piston=petal_piston, \
            atm_screen=atm_screen, ncpa_screen=ncpa_screen, **conf)
    conf['RAVC_misalign'] = misalign
    # define flags for special cases
    RAVC = True if conf['mode'] in ['RAVC'] else False
    APP = True if conf['mode'] in ['APP'] and conf['onaxis'] is True else False 
    # pupil-plane apodization
    wf, apo_amp, apo_phase = apodization(wf, conf, RAVC=RAVC)
    # select focal-plane mask, and Lyot-stop
    if conf['mode'] in ['ELT', 'APP']:
        conf['LS_params'] = [1., -0.3, 0.] # no Lyot stop
    elif conf['mode'] in ['CLC']:
        conf['LS_params'] = [0.8, 0.1, 1.1]
        if conf['onaxis'] == True:
            lyot(wf, conf)
    elif conf['mode'] in ['CVC', 'RAVC']:
        conf['LS_params'] = [0.98, 0.03, 1.1]
        if conf['onaxis'] == True:
            vortex(wf, conf)
    # lyot-stop
    wf, LS_amp, LS_phase = lyotstop(wf, conf, RAVC=RAVC, APP=APP)
    # get science image
    psf = detector(wf, conf)
    # release memory: global variable set to None
    wf = None
    
    if full_output is True:
        return psf, LS_amp, apo_amp
    else:
        return psf

""" load all wavefront aberrations """

atm_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))[:3]
#atm_cube = fits.getdata('/mnt/disk4tb/METIS/COMPASS_243x243/cube_COMPASS_20181008_3600s_100ms.fits')[::3,:]
nframes = atm_cube.shape[0]
print('\nNumber of frames = %s.'%nframes)
ncpa_cube = fits.getdata(os.path.join(conf['input_dir'], conf['ncpa_screen_file']))
ncpa_scaling = 10*0.0089
ncpa_cube = np.array([x*ncpa_cube for x in np.linspace(-ncpa_scaling,ncpa_scaling,nframes)])
point_drift = (np.linspace(-0.2, 0.2, nframes)*u.mas).to('rad').value/(3.8e-6/37)
point_jitter = (np.random.normal(0, 1, nframes)*u.mas).to('rad').value/(3.8e-6/37)
pupil_drift = [[x,0,0,0,0,0] for x in np.linspace(-0.01,0.01,nframes)]
np.random.seed(345678)
piston_drift_ptv = 0.01 #µm 
piston_drift_freq = 1 #1/T_adi
piston_drift = np.random.normal(piston_drift_ptv/2,piston_drift_ptv/2/10,6) \
        *np.sin(np.array([range(nframes)]).T/nframes \
        *np.random.normal(2*np.pi*piston_drift_freq,2*np.pi*piston_drift_freq/10,6) \
        + np.random.uniform(0,2*np.pi,6))
cases = [[0, atm_cube, ncpa_cube, None        , None         , None                    , 0, 0],
         [1, atm_cube, None     , None        , None         , None                    , 0, 7],
         [2, atm_cube, None     , None        , None         , point_drift             , 2, 0],
         [3, atm_cube, None     , None        , None         , point_jitter            , 2, 0],
         [4, atm_cube, None     , None        , pupil_drift  , None                    , 0, 0],
         [5, atm_cube, None     , piston_drift, None         , None                    , 0, 0],
         [6, atm_cube, ncpa_cube, piston_drift, pupil_drift  , point_drift+point_jitter, 2, 7]]
#cases = [cases[i] for i in [2,3,4,5,6]]

""" Start looping on the different cases, bands, modes """

for case in cases:
    atm_screens = case[1] if np.any(case[1]) else [None]*nframes
    ncpa_screens = case[2] if np.any(case[2]) else [None]*nframes
    petal_pistons = case[3] if np.any(case[3]) else [None]*nframes
    misaligns = case[4] if np.any(case[4]) else [None]*nframes
    zernikes = case[5] if np.any(case[5]) else [None]*nframes
    conf['zern_inds'] = case[6]
    conf['N_mis_segments'] = case[7]

    for band in conf['bands']:
        conf['band'] = band
        conf['lam'] = band_specs[band]['lam']
        conf['pscale'] = band_specs[band]['pscale']
    
        # compute beam ratio, pupil size, and create the entrance pupil
        wf_start, pup_amp, pup_phase = pupil(conf)
    
        for i, mode in enumerate(band_specs[band]['modes']):
            conf['mode'] = mode
        
            # starting time
            t0 = time.time()
        
            # propagate the frames, using multiple cores if possible 
            if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
                if conf['cpucount'] == None:
                    conf['cpucount'] = mpro.cpu_count() - 1
                print('%s: %s band, %s mode, using %s cores.'\
                        %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                        band, mode, conf['cpucount']))
                p = mpro.Pool(conf['cpucount'])
                func = partial(propagate, wf_start, conf, False)
                psfs = np.array(p.starmap(func, zip(atm_screens, ncpa_screens, petal_pistons, misaligns, zernikes)))
            else:
                print('%s: %s band, %s mode, using %s core.'\
                        %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                        band, mode, 1))
                psfs = np.zeros((nframes, conf['ndet'], conf['ndet']))
                for atm_screen, ncpa_screen, petal_piston, misalign, zernike in \
                        zip(atm_screens, ncpa_screens, petal_pistons, misaligns, zernikes):
                    psf, LS_amp, apo_amp = propagate(wf_start, conf, True, \
                            atm_screen, ncpa_screen, petal_piston, misalign, zernike)
                    psfs[i,:,:] = psf
            
            # if only one frame, make dim = 2
            if nframes == 1:
                psfs = psfs[0]
            
            """ Write to .fits """
            
            # file name
            conf['prefix'] = '%s_'%case[0]
            on_off = {True:'onaxis', False:'offaxis'}[conf['onaxis']]
            filename = '%s%s'%(conf['prefix'], on_off)+'_%s_'+'%s_%s'%(band, mode)
            # save PSF (or PSFs cube) and/or Lyot-Stop
            fits.writeto(os.path.join(conf['output_dir'], filename%'PSF') \
                    + '.fits', np.float32(psfs), overwrite=True)
            if False:
                fits.writeto(os.path.join(conf['output_dir'], filename%'APO') \
                        + '.fits', apo_amp, overwrite=True)
                fits.writeto(os.path.join(conf['output_dir'], filename%'PUP') \
                        + '.fits', pup_amp, overwrite=True)
                fits.writeto(os.path.join(conf['output_dir'], filename%'LS') \
                        + '.fits', LS_amp, overwrite=True)
            # release memory
            psfs = None
            
            
            """ Save figures to .png """
            
            if False:
                plt.figure()
                #plt.imshow(psf**0.05, origin='lower')
                plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
                plt.colorbar()
                plt.show(block=False)
                plt.savefig(os.path.join(conf['output_dir'], filename%'PSF') \
                        + '.png', dpi=300, transparent=True)
            if False:
                plt.figure()
                plt.imshow(pup_amp[50:-50,50:-50], origin='lower', cmap='gray', vmin=0, vmax=1)
                #plt.imshow(pup_amp[50:-50,50:-50], origin='lower')
                plt.colorbar()
                plt.show(block=False)
                plt.savefig(os.path.join(conf['output_dir'], filename%'LS') \
                        + '.png', dpi=300, transparent=True)
            
            # print elapsed time
            print('      Elapsed %.3f seconds.'%(time.time() - t0))

# send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))
