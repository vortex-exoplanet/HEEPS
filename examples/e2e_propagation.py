#
# End-to-end propagation with HEEPS
#
# AIM:
#   script file to configure and launch one or several simulations.
#   Each simulation generates a post-coronagraphic PSF cube saved in the output
#   directory.
#
# TODO:
#    - create a MySimulation class with the contributors configuration, the propagate function, etc. all in one
#    - improve diagnostics/logging saving
#    - move band_specs to config.py
#
# HISTORY:
#  - 2020-03-30: implementing configuration file and contributor saving. (ULiege)
#  - 2020-03-27: refurbished by GOX  based on CD version (ULiege)

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

from util_e2e_propagation import loopOverCasesBandsModes, saveHeepsConfiguration

#----------------------------------#
# * General script variables *
load_effects = True  # jitter, ncpa, etc.
# activate saving of this script file, the conf dic, and the contributor script
logging = True
# output dir where the post-coronagraphic PSF will be save as well as the logging info
output_dir = '/mnt/disk4tb/METIS/heeps_analysis/output_files_gox/20200330/'

# create directory if it does not exists
if not(os.path.isdir(output_dir)):
    print('Creating output directory %s'.format(output_dir))
    os.mkdir(output_dir)

if logging:
    print('Copying script file to output directory')
    today = datetime.now()
    time_tag = today.strftime('%Y%m%d_%H%M%S')
    copyfile(__file__,
             output_dir + time_tag + '_' +\
             os.path.basename(__file__))

#----------------------------------#
# * Changing HEEPS configuration *
conf['output_dir'] = output_dir
conf['spiders_width'] = 0.5
conf['LS_params'] = [0.98, 0.03, 0]
conf['bands'] = ['L'] #, 'M', 'N1', 'N2', 'N1a', 'N2a']
conf['CLC_diam'] = 4
conf['onaxis'] = False
conf['static_ncpa'] = False
conf['send_to'] = 'gorban@uliege.be'
conf['send_message'] = 'End-to-end simulation finished OK.'
conf['get_amp'] = False
conf['get_phase'] = False

# specify CPU count: 1=single processing, None=max-1
conf['cpucount'] = 28

# * Band specifications
# TODO: move band_specs to config.py
# modes, should be a separate variable, which includes standards Lyot stop parameter?

band_specs = {'L': {'lam': 3.812e-6,# 3.8204e-6,#
                  'pscale': 5.21,
                   'modes': ['RAVC']},#['ELT', 'RAVC', 'CVC', 'APP', 'CLC']},
             }

# * Sim configuration
conf['hfov'] = 0.666   # 1.3# #half field of view in arcsec
conf['gridsize'] = 1024
# conf['pupil_file'] = 'ELT_pupil_downsampled_Lp.fits'
# conf['pupil_file'] = 'rot90_M_ELT_binary_pupil.fits'
# conf['pupil_file'] = 'ELT_37_0.3_1025.fits'
conf['pupil_file'] = ''  # no spiders; no segments

# Saving configuration file
# Do not change 'conf' after
if logging:
    print('Copying the config OrderedDic to .npy')
    suffix = ''
    conf_fname = conf['output_dir'] + time_tag + '_heeps_config' + suffix + '.npy'
    saveHeepsConfiguration(conf_fname, conf)

#----------------------------------#
#  ALL CONTRIBUTORS

# * Load all wavefront aberrations
data_dir = '/mnt/disk4tb/METIS/COMPASS_GOX/'
phase_screen_file = data_dir + \
                'cube_COMPASS_20190430_600s_100ms_lag_2_g0.4_noPiston_noTT.fits'
atm_cube = fits.getdata(phase_screen_file)

#nframes = atm_cube.shape[0]
nframes = 1   # for ELT e.g.
print('\nNumber of frames = %s.'%nframes)

if load_effects:
    from set_contributors  import *
    if logging:
        print('Copying the script file that configure all the effects')
        copyfile(os.path.dirname(__file__) + '/' + 'set_contributors.py',
                 conf['output_dir'] + time_tag + '_' + 'set_contributors.py')

else:
    ncpa_piston_ALL = None
    pupil_drift = None
    point_ALL = None

#----------------------------------------#
# CASES :
# Entries gives
#   - [0] id number
#   - [1] Atmospheric turbulence cube
#          - units: ...
#   - [2] NCPA cube
#          - units: nm rms
#          - STA, QLSF, QHSF, DYN, ALL, +piston_all
#   - [3] Piston
#          - units: nm rms
#          - QLSF, QHSF, DYN, do not used if using ncpa_piston_all
#   - [4] Pupil drift
#          - units: fraction of pupil (0.01 means 1%)
#   - [5] Zernike coefficients
#          - units: mas or in lambda / D, or nm rms ??
#          - QSTA, DYN, ALL
#   - [6] Zernike mode indices
#          - [2,3] for tip and tilt
#   - [7] Numbers of missing segments
#
#  To inject pointing jitter, use [6] and [7]

cases = [[0,  None,     None       , None        , None         , None         , 0, 0],
         [1,  atm_cube, None       , None        , None         , None         , 0, 0],
         [81,  atm_cube, ncpa_piston_ALL, None, pupil_drift    , point_ALL    , [2,3], 0],
         [None]]

cases = [case for case in cases if case[0] in [0]]
print('All cases created: %s'%[cases[0] for case in cases])


#----------------------------------------#
# >> Run simulations !
# or Create your own loop here

loopOverCasesBandsModes(conf, nframes, cases, band_specs)


#----------------------------------------#
# send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))
