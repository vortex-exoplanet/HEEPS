#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
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

""" user inputs """

conf['bands'] = ['L']# , 'M', 'N1', 'N2']
conf['CLC_diam'] = 4
conf['onaxis'] = True
conf['static_ncpa'] = False
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_message'] = 'End-to-end simulation finished OK.'
conf['get_amp'] = False
conf['get_phase'] = False
conf['APP_phase_file'] = 'rot90_METIS_vAPP_PDR_L.fits'# 'rot90_METIS_vAPP_PDR_M.fits'#

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = None

# band specifications
band_specs = {'L': {'lam': 3.812e-6,# 3.8204e-6,#
                  'pscale': 5.21,
                   'modes': ['RAVC']},#['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
              'M': {'lam': 4.8091E-06,# 4.7970e-6,#
                  'pscale': 5.21,
                   'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
              'N1': {'lam': 8.6313e-6, #8.6016e-6,
                  'pscale': 6.79,
                   'modes': ['ELT', 'CVC', 'CLC']},
              'N2': {'lam': 11.287e-6, #11.236e-6,
                  'pscale': 6.79,
                   'modes': ['ELT', 'CVC', 'CLC']}}

""" Create a function to propagate one single wavefront """

def propagate(wf_start, conf, full_output, atm_screen, ncpa_screen, 
            piston_screen, misalign, zernike):
    # keep a copy of the input pupil wavefront
    wf = deepcopy(wf_start)
    # apply wavefront aberrations
    wavefront_aberrations(wf, zernike=zernike, piston_screen=piston_screen, \
            atm_screen=atm_screen, ncpa_screen=ncpa_screen, **conf) #resized=True, nans=False, **conf)
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
    conf['ndet'] = int(np.ceil(2*conf['hfov']*1000/conf['pscale']))    
    if conf['ndet'] % 2:
        conf['ndet'] += 1
    psf = detector(wf, conf)
    # release memory: global variable set to None
    wf = None
    
    if full_output is True:
        return psf, LS_amp, apo_amp
    else:
        return psf

""" load all wavefront aberrations """

conf['hfov'] = 0.666# 1.3#
conf['gridsize'] = 1024
conf['pupil_file'] = 'ELT_37_0.3_1025.fits'
#conf['atm_screen_file'] = 'L_band_253/cube_COMPASS_20181008_3600s_300ms.fits'
atm_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))[:3]
nframes = atm_cube.shape[0]# 1#
print('\nNumber of frames = %s.'%nframes)

if False:
    # load ncpa maps (spatial frequencies)    
    ncpa_allSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_allSF_253.fits'))
    ncpa_LSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_LSF_10cpp_253.fits'))
    ncpa_HSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_HSF_10cpp_253.fits'))
    mean_allSF = np.nanmean(ncpa_allSF)
    mean_LSF = np.nanmean(ncpa_LSF)
    mean_HSF = np.nanmean(ncpa_HSF)

    # load time series (temporal frequencies)
    LTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=832404.fits'))
    LTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=523364.fits'))
    LTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=409566.fits'))
    LTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=224788.fits'))
    HTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=832404.fits'))
    HTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=523364.fits'))
    HTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=409566.fits'))
    HTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=224788.fits'))

    # create ncpa cubes
    #ncpa_STA = np.array([ncpa_HSF]*nframes)
    ncpa_STA = np.array([ncpa_HSF]*12000)
    ncpa_QLSF = np.array([x*ncpa_LSF for x in LTF1])
    ncpa_QHSF = np.array([x*ncpa_HSF for x in LTF2])
    ncpa_DYN = np.array([x*ncpa_allSF for x in HTF1])
    ncpa_ALL = ncpa_STA*35.9 + ncpa_QLSF*20 + ncpa_QHSF*20 + ncpa_DYN*40

    # petal piston
    piston_QLSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_LSF_1rms_seed=123456.fits'))
    piston_QHSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_HSF_1rms_seed=234567.fits'))
    piston_DYN = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_HTF_allSF_1rms_seed=345678.fits'))

    # ncpa + piston:
    ncpa_piston_QLSF = (ncpa_QLSF + piston_QLSF)/np.sqrt(2)
    ncpa_piston_QHSF = (ncpa_QHSF + piston_QHSF)/np.sqrt(2)
    ncpa_piston_DYN = (ncpa_DYN + piston_DYN)/np.sqrt(2)
    # norm (almost 1)
    ncpa_piston_QLSF /= np.mean([np.nanstd(x) for x in ncpa_piston_QLSF])
    ncpa_piston_QHSF /= np.mean([np.nanstd(x) for x in ncpa_piston_QHSF])
    ncpa_piston_DYN /= np.mean([np.nanstd(x) for x in ncpa_piston_DYN])
    # total
    ncpa_piston_ALL = ncpa_STA*35.9 + ncpa_piston_QLSF*20 + ncpa_piston_QHSF*20 + ncpa_piston_DYN*40

    # apodizer drift
    drift_ptv = 0.01 # 1% ptv
    pupil_drift = np.array([[x,0,0,0,0,0] for x in np.linspace(-drift_ptv/2, drift_ptv/2, nframes)])

    # create pointing errors (zernikes [2,3])
    #point_drift_x = np.linspace(-0.2, 0.2, nframes)
    #point_jit_x = (np.random.normal(0, 2, nframes)*u.mas).to('rad').value/(3.8e-6/37)
    point_QSTA = np.array([LTF3, LTF4]).T/np.sqrt(2) # in xy
    point_DYN = np.array([HTF3, HTF4]).T/np.sqrt(2) # in xy
    point_ALL = point_QSTA*0.4 + point_DYN*2      # factors in mas  


# cases
cases = [[0,  None,     None       , None        , None         , None         , 0, 0],
         [1,  atm_cube, None       , None        , None         , None         , 0, 0],

         #[11,  atm_cube, ncpa_piston_ALL, None, pupil_drift    , point_ALL    , [2,3], 0],
         #[22,  atm_cube, ncpa_piston_ALL, None, pupil_drift*2  , point_ALL    , [2,3], 0],
         #[33,  atm_cube, ncpa_piston_ALL, None, pupil_drift*3  , point_ALL    , [2,3], 0],

         [None]]

cases = [case for case in cases if case[0] in [0,1]]
#cases = [case for case in cases if case[0] in [11,12,13]]
print('Case numbers = %s'%[case[0] for case in cases])

""" Start looping on the different cases, bands, modes """

for case in cases:
    atm_screens = case[1] if np.any(case[1]) else [None]*nframes
    ncpa_screens = case[2] if np.any(case[2]) else [None]*nframes
    piston_screens = case[3] if np.any(case[3]) else [None]*nframes
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
        
        for mode in band_specs[band]['modes']:
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
                psfs = np.array(p.starmap(func, zip(atm_screens, ncpa_screens, \
                        piston_screens, misaligns, zernikes)))
                p.close()
                p.join()
            else:
                print('%s: %s band, %s mode, using %s core.'\
                        %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                        band, mode, 1))
                psfs = np.zeros((nframes, conf['ndet'], conf['ndet']))
                for i, (atm_screen, ncpa_screen, piston_screen, misalign, \
                        zernike) in enumerate(zip(atm_screens, ncpa_screens, \
                        piston_screens, misaligns, zernikes)):
                    psf, LS_amp, apo_amp = propagate(wf_start, conf, True, \
                            atm_screen, ncpa_screen, piston_screen, misalign, zernike)
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
