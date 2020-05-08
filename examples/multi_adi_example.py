#import matplotlib; matplotlib.use('agg') # to run on a headless server
from heeps.config import conf
from heeps.contrast.adi import adi
import multiprocessing as mpro
from functools import partial
from sys import platform
import time
import os
from heeps.config import conf

# email when finished
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_message'] = 'ADI simulation finished OK.'

# list of cases
cases = [{'band':'N1', 'mode':'CVC','pscale':6.79, 'mag': 6,  'add_bckg':False, 'folder':'scao_only'},
         {'band':'N1', 'mode':'CVC','pscale':6.79, 'mag': 6,  'add_bckg':False, 'folder':'all_effects'}]

# number of cases
ncases = len(cases)
# number of CPUs to use: one case per CPU
conf['cpucount'] = ncases

# define a wrapper function
def multi_cc(verbose, case, **conf):
    cube_duration = 3600        # used for 'savename'
    cube_samp = 300             # used for 'savename'
    band = case['band']
    mode = case['mode']
    pscale = case['pscale']
    add_bckg = case['add_bckg']
    mag = case['mag']
    # optional sub-folder for onaxis PSFs
    folder = str(case.get('folder',''))
    tagname = '_%s'%folder
#    path_onaxis = 'cube_COMPASS_20181008_3600s_300ms_12000x256x256/app/%s/'%folder
#    path_onaxis = 'cube_COMPASS_201910_Gilles/%s/'%folder
    path_onaxis = 'output_files/%s'%folder
#    path_offaxis = 'offaxis_2048_1.3'    
#    path_offaxis = 'offaxis_256x256'    
    path_offaxis = path_onaxis
    path_output = path_onaxis
    # call adi function
    adi(path_offaxis=path_offaxis, path_onaxis=path_onaxis, \
            path_output=path_output, verbose=verbose, \
            mode=mode, band=band, cube_duration=cube_duration, \
            cube_samp=cube_samp, pscale=pscale, \
            add_bckg=add_bckg, mag=mag, prefix=folder, tagname=tagname,
            **conf)
    # print stuff
    print('%s: Finished case %s.'%(time.strftime("%Y-%m-%d %H:%M:%S", \
            time.localtime()), folder))

# start ADI
t0 = time.time()
print('\n%s: ADI simulation starts, %s cases.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ncases))
if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
    print('Multiprocessing using %s cores'%conf['cpucount'])
    verbose = False
    func = partial(multi_cc, verbose, **conf)
    p = mpro.Pool(conf['cpucount'])
    p.map(func, cases)
    p.close()
    p.join()
else:
    for case in cases:
        verbose = True
        multi_cc(verbose, case)

# print elapsed time
print('      Elapsed %.3f seconds.'%(time.time() - t0))
print(time.strftime("%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))

# send email when simulation finished
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))
