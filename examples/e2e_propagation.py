import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import fits 
from heeps.config import conf
from heeps.pupil import pupil
import heeps.util.img_processing as impro
import heeps.util.propagate as prop
import os.path
import time
import multiprocessing as mpro
from functools import partial
from sys import platform
import os

""" user inputs """

conf['bands'] = ['N1']#, 'M', 'N1', 'N2', 'N1a', 'N2a']
conf['CLC_diam'] = 4
conf['onaxis'] = True
conf['static_ncpa'] = False
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_message'] = 'End-to-end simulation finished OK.'
conf['get_amp'] = False
conf['get_phase'] = False
conf['APP_phase_file'] = 'rot90_METIS_vAPP_PDR_L.fits'# 'rot90_METIS_vAPP_PDR_M.fits'#
conf['hfov'] = 0.666#1.3# 
conf['nframes'] = 50

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = None

# band specifications
band_specs = {'L': {'lam': 3.8113E-06,
                  'pscale': 5.47,
                   'modes': ['RAVC']},#['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
              'M': {'lam': 4.8153E-06,
                  'pscale': 5.47,
                   'modes': ['RAVC']},
              'N1': {'lam': 8.6904e-6,
                  'pscale': 6.79,
                   'modes': ['CVC']},
              'N2': {'lam': 11.242e-6,
                  'pscale': 6.79,
                   'modes': ['CVC']},
              'N1a':{'lam': 8.6658e-6, #8.6469e-6,
                  'pscale': 10.78,
                   'modes': ['ELT', 'CVC', 'CLC']},
              'N2a':{'lam': 11.219e-6, #11.187e-6,
                  'pscale': 10.78,
                   'modes': ['CVC']}}#'ELT', 'CVC', 'CLC']}}

""" load all wavefront aberrations """

conf['atm_screen_file'] = 'cube_COMPASS_20181008_3600s_300ms.fits'
atm_cube = fits.getdata(os.path.join(conf['input_dir'], conf['atm_screen_file']))
# number of frames
if conf['nframes'] is None:
    conf['nframes'] = len(atm_cube)
else:
    atm_cube = atm_cube[:conf['nframes']]
print('nframes = %s'%conf['nframes'])

if True:
    # load ncpa maps (spatial frequencies)    
    ncpa_allSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_allSF_253.fits'))
    ncpa_LSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_LSF_10cpp_253.fits'))
    ncpa_HSF = fits.getdata(os.path.join(conf['input_dir'], 'ncpa_HSF_10cpp_253.fits'))
    mean_allSF = np.nanmean(ncpa_allSF)
    mean_LSF = np.nanmean(ncpa_LSF)
    mean_HSF = np.nanmean(ncpa_HSF)

    # load time series (temporal frequencies)
    LTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=832404.fits'))[:conf['nframes']]
    LTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=523364.fits'))[:conf['nframes']]
    LTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=409566.fits'))[:conf['nframes']]
    LTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_LTF_0-0.01Hz_12000x1rms_seed=224788.fits'))[:conf['nframes']]
    HTF1 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=832404.fits'))[:conf['nframes']]
    HTF2 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=523364.fits'))[:conf['nframes']]
    HTF3 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=409566.fits'))[:conf['nframes']]
    HTF4 = fits.getdata(os.path.join(conf['input_dir'], \
            'time_series_HTF_0.01-1Hz_12000x1rms_seed=224788.fits'))[:conf['nframes']]

    # create ncpa cubes
    ncpa_STA = np.array([ncpa_HSF]*conf['nframes'])
    ncpa_QLSF = np.array([x*ncpa_LSF for x in LTF1])
    ncpa_QHSF = np.array([x*ncpa_HSF for x in LTF2])
    ncpa_DYN = np.array([x*ncpa_allSF for x in HTF1])
    ncpa_ALL = ncpa_STA*35.9 + ncpa_QLSF*20 + ncpa_QHSF*20 + ncpa_DYN*40

    # petal piston
    piston_QLSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_LSF_1rms_seed=123456.fits'))[:conf['nframes']]
    piston_QHSF = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_LTF_HSF_1rms_seed=234567.fits'))[:conf['nframes']]
    piston_DYN = fits.getdata(os.path.join(conf['input_dir'], \
            'cube_petal_piston_HTF_allSF_1rms_seed=345678.fits'))[:conf['nframes']]

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
    pupil_drift = np.array([[x,0,0,0,0,0] for x in np.linspace(-drift_ptv/2, drift_ptv/2, conf['nframes'])])

    # create pointing errors (zernikes [2,3])
    point_QSTA = np.array([LTF3, LTF4]).T/np.sqrt(2) # in xy
    point_DYN = np.array([HTF3, HTF4]).T/np.sqrt(2) # in xy
    point_ALL = point_QSTA*0.4 + point_DYN*2      # factors in mas  


# cases
cases = [['perfect',  None,     None       , None        , None         , None         , 0, 0],
         ['scao_only',  atm_cube, None       , None        , None         , None         , 0, 0],
         ['all_effects',  atm_cube, ncpa_piston_ALL, None, pupil_drift*3  , point_ALL    , [2,3], 0],
         ['mis_seg',  atm_cube, ncpa_piston_ALL, None, pupil_drift*3  , point_ALL    , [2,3], 7],

         [None]]
cases = [cases[i] for i in [1,2]]
print('selected cases: %s'%[case[0] for case in cases])

""" Start looping on the different cases, bands, modes """

for case in cases:
    # create case folders
    case_folder = os.path.join(conf['output_dir'], '%s'%case[0])
    os.makedirs(case_folder, exist_ok=True)
    print("\nCreating case folder: '%s'"%case[0])

    atm_screens = case[1] if np.any(case[1]) else [None]*conf['nframes']
    ncpa_screens = case[2] if np.any(case[2]) else [None]*conf['nframes']
    piston_screens = case[3] if np.any(case[3]) else [None]*conf['nframes']
    misaligns = case[4] if np.any(case[4]) else [None]*conf['nframes']
    zernikes = case[5] if np.any(case[5]) else [None]*conf['nframes']
    conf['zern_inds'] = case[6]
    conf['N_mis_segments'] = case[7]
    
    for band in conf['bands']:
        conf['band'] = band
        conf['lam'] = band_specs[band]['lam']
        conf['pscale'] = band_specs[band]['pscale']

        # compute beam ratio, pupil size, and create the entrance pupil
        conf['full_output'] = True
        wf_start, pup_amp, pup_phase = pupil(conf)
        if np.any(pup_amp):
            fits.writeto(os.path.join(case_folder, '%spupil_%s.fits'\
                %(conf['prefix'], band)), np.float32(pup_amp), overwrite=True)
        print('   band=%s, npupil=%s'%(band, conf['npupil']))
            
        # create ELT off-axis PSF
        conf['mode'] = 'ELT'
        conf['onaxis'] = False
        psf, _, _ = prop.prop_one(wf_start, conf, \
                            atm_screens[0], ncpa_screens[0], piston_screens[0], misaligns[0], zernikes[0])
        fits.writeto(os.path.join(case_folder, \
            '%soffaxis_PSF_%s_ELT.fits'%(conf['prefix'], band)), np.float32(psf), overwrite=True)

        for mode in band_specs[band]['modes']:
            # create off-axis PSF for selected mode
            conf['full_output'] = True
            conf['mode'] = mode
            conf['onaxis'] = False
            psf, apo_amp, LS_amp = prop.prop_one(wf_start, conf, \
                            atm_screens[0], ncpa_screens[0], piston_screens[0], misaligns[0], zernikes[0])
            fits.writeto(os.path.join(case_folder, \
                '%soffaxis_PSF_%s_%s.fits'%(conf['prefix'], band, mode)), np.float32(psf), overwrite=True)
            if np.any(apo_amp):
                fits.writeto(os.path.join(case_folder, '%spupil_%s_%s_apodizer.fits'\
                    %(conf['prefix'], band, mode)), np.float32(apo_amp), overwrite=True)
            if np.any(LS_amp):
                fits.writeto(os.path.join(case_folder, '%spupil_%s_%s_LS_offaxis.fits'\
                    %(conf['prefix'], band, mode)), np.float32(LS_amp), overwrite=True)

            # propagate 1 on-axis PSF for selected mode
            conf['onaxis'] = True
            t0 = time.time()
            psf, apo_amp, LS_amp = prop.prop_one(wf_start, conf, \
                            atm_screens[0], ncpa_screens[0], piston_screens[0], misaligns[0], zernikes[0])
            print('      mode=%s, estimated duration: %.3f seconds.'%(mode, (time.time() - t0)*conf['nframes']))
            if np.any(LS_amp):
                fits.writeto(os.path.join(case_folder, '%spupil_%s_%s_LS_onaxis.fits'\
                    %(conf['prefix'], band, mode)), np.float32(LS_amp), overwrite=True)
    
            # propagate onaxis cube, using multiple cores if possible 
            conf['full_output'] = False            
            t0 = time.time()
            if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
                if conf['cpucount'] == None:
                    conf['cpucount'] = mpro.cpu_count() - 1
                print('      %s: %s band, %s mode, using %s cores.'\
                        %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                        band, mode, conf['cpucount']))
                p = mpro.Pool(conf['cpucount'])
                func = partial(prop.prop_one, wf_start, conf)
                psfs = np.array(p.starmap(func, zip(atm_screens, ncpa_screens, \
                        piston_screens, misaligns, zernikes)))
                p.close()
                p.join()
            else:
                print('      %s: %s band, %s mode, using %s core.'\
                        %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                        band, mode, 1))
                psfs = np.zeros((conf['nframes'], conf['ndet'], conf['ndet']))
                for i, (atm_screen, ncpa_screen, piston_screen, misalign, \
                        zernike) in enumerate(zip(atm_screens, ncpa_screens, \
                        piston_screens, misaligns, zernikes)):
                    psf, LS_amp, apo_amp = prop.prop_one(wf_start, conf, \
                            atm_screen, ncpa_screen, piston_screen, misalign, zernike)
                    psfs[i,:,:] = psf                
            # if only one frame, make dim = 2
            if conf['nframes'] == 1:
                psfs = psfs[0]
            # save onaxis PSFs
            fits.writeto(os.path.join(case_folder, '%sonaxis_PSF_%s_%s.fits'\
                %(conf['prefix'], band, mode)), np.float32(psfs), overwrite=True)
            # release memory
            psfs = None            
            # print elapsed time
            print('      elapsed %.3f seconds.'%(time.time() - t0))

# send email when simulation finished
print(time.strftime("\n%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))

from heeps.contrast import contrast_raw
import multi_adi_example


%pylab
import astropy.units as u
#import examples.multi_adi_example

#lam = 8.6904e-6
#diam = 36.3508
lamD = conf['lam']/conf['diam']*u.rad.to('arcsec')
raw1 = fits.getdata('scao_only/cc_raw_N1_CVC_scao_only.fits')
raw2 = fits.getdata('all_effects/cc_raw_N1_CVC_all_effects.fits')
adi1= fits.getdata('scao_only/cc_compass3600s_samp300ms_dec-5_N1mag6_bckg0_CVC_scao_only.fits')
adi2= fits.getdata('all_effects/cc_compass3600s_samp300ms_dec-5_N1mag6_bckg0_CVC_all_effects.fits')

figure(figsize=(12,4))
plot(raw2[0,:],raw2[1,:], label='raw, all effects')
plot(raw1[0,:],raw1[1,:], label='raw, scao only')
plot(adi2[:,4]/lamD,adi2[:,1],label='adi, all effects')
plot(adi1[:,4]/lamD,adi1[:,1],label='adi, scao only')

loglog()
legend()
grid(True)
xlim(1,30)
ylim(1e-8,1e-1)

savefig('test_50.png')