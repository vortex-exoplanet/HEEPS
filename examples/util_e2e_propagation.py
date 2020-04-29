import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import fits
from heeps.config import conf
from heeps.pupil import pupil
from heeps.aberrations import wavefront_aberrations
from heeps.coronagraphs import apodization, vortex, lyotstop, lyot
from heeps.detector import detector
import heeps.util.img_processing as impro
import os.path
from copy import deepcopy
import time
import multiprocessing as mpro
from functools import partial
from sys import platform
import os
from shutil import copyfile
from datetime import datetime

def saveHeepsConfiguration(confFilename, conf):
    assert confFilename.endswith('.npy')
    np.save(confFilename, conf)

def loadHeepsConfigurationFromFile(confFilename):
    assert confFilename.endswith('.npy')
    return np.load(confFileName).item()

def _checkIfLsIsStandard(ls_std, ls_conf, mode=None):
    if not(ls_std == ls_conf):
        print('Warning, Lyot Stop parameters not standard for mode %s'%mode)

def propagate(wf_start, conf, full_output, atm_screen, ncpa_screen,
            piston_screen, misalign, zernike):
    """ Create a function to propagate one single wavefront """

    # * keep a copy of the input pupil wavefront
    wf = deepcopy(wf_start)
    # * apply wavefront aberrations
    wavefront_aberrations(wf, zernike=zernike, piston_screen=piston_screen, \
            atm_screen=atm_screen, ncpa_screen=ncpa_screen, **conf)
#            atm_screen=atm_screen, ncpa_screen=ncpa_screen, resized=True, nans=False, **conf)
    conf['RAVC_misalign'] = misalign
    # * define flags for special cases
    RAVC = True if conf['mode'] in ['RAVC'] else False
    APP = True if conf['mode'] in ['APP'] and conf['onaxis'] is True else False
    # * pupil-plane apodization
    wf, apo_amp, apo_phase = apodization(wf, conf, RAVC=RAVC)
    # * select focal-plane mask, and Lyot-stop
    if conf['mode'] in ['ELT', 'APP']:
        #conf['LS_params'] = [1., -0.3, 0.] # no Lyot stop

        # TODO: explain why '-0.3' for dR_in(%) here
        ls_standard = [1., -0.3, 0.] # no Lyot stop
        _checkIfLsIsStandard(ls_standard, conf['LS_params'], mode=conf['mode'])
    elif conf['mode'] in ['CLC']:
        #conf['LS_params'] = [0.8, 0.1, 1.1]
        ls_standard = [0.8, 0.1, 1.1]
        _checkIfLsIsStandard(ls_standard, conf['LS_params'], mode=conf['mode'])
        if conf['onaxis'] == True:
            lyot(wf, conf)
    elif conf['mode'] in ['CVC', 'RAVC']:
        #conf['LS_params'] = [0.98, 0.03, 1.1]
        ls_standard = [0.98, 0.03, 1.1]
        _checkIfLsIsStandard(ls_standard, conf['LS_params'], mode=conf['mode'])
        if conf['onaxis'] == True:
            vortex(wf, conf)
    # * lyot-stop
    wf, LS_amp, LS_phase = lyotstop(wf, conf, RAVC=RAVC, APP=APP)
    # * get science image
    conf['ndet'] = int(np.ceil(2*conf['hfov']*1000/conf['pscale']))
    if conf['ndet'] %2 :
        conf['ndet'] += 1
    psf = detector(wf, conf)
    # release memory: global variable set to None
    wf = None

    if full_output is True:
        return psf, LS_amp, apo_amp
    else:
        return psf


def loopOverCasesBandsModes(conf, nframes, cases, band_specs,
                            suffix=None, savefigs=False, overwrite=False):
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
                filename = '%s%s'%(conf['prefix'], on_off)+'_%s_'+'%s_%s%s'%(band, mode, suffix)

                filenamePSF = os.path.join(conf['output_dir'], filename%'PSF') + '.fits'

                if os.path.isfile(filenamePSF) and not(overwrite):
                    print('Warning, file %s does already exist, appending time tag'%filenamePSF)
                    today = datetime.now()
                    time_tag = today.strftime('%Y%m%d_%H%M%S')
                    filenamePSF = conf['output_dir'] + os.path.basename(filenamePSF) + \
                        '_' + time_tag + '.fits'

                # save PSF (or PSFs cube) and/or Lyot-Stop
                fits.writeto(filenamePSF, np.float32(psfs), overwrite=overwrite)
#                 if False:
#                     fits.writeto(os.path.join(conf['output_dir'], filename%'APO') \
#                             + '.fits', apo_amp, overwrite=False)
#                     fits.writeto(os.path.join(conf['output_dir'], filename%'PUP') \
#                             + '.fits', pup_amp, overwrite=False)
#                     fits.writeto(os.path.join(conf['output_dir'], filename%'LS') \
#                             + '.fits', LS_amp, overwrite=False)

                # release memory
                psfs = None


                """ Save figures to .png """

                if savefigs:
                    plt.figure()
                    #plt.imshow(psf**0.05, origin='lower')
                    plt.imshow(np.log10(psf/1482.22), origin='lower') # 1482.22 is peak in ELT mode
                    plt.colorbar()
                    plt.show(block=False)
                    plt.savefig(os.path.join(conf['output_dir'], filename%'PSF') \
                            + '.png', dpi=300, transparent=True)
                if savefigs:
                    plt.figure()
                    plt.imshow(pup_amp[50:-50,50:-50], origin='lower', cmap='gray', vmin=0, vmax=1)
                    #plt.imshow(pup_amp[50:-50,50:-50], origin='lower')
                    plt.colorbar()
                    plt.show(block=False)
                    plt.savefig(os.path.join(conf['output_dir'], filename%'LS') \
                            + '.png', dpi=300, transparent=True)

                # print elapsed time
                print('      Elapsed %.3f seconds.'%(time.time() - t0))
