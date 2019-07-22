#import matplotlib; matplotlib.use('agg') # to run on a headless server
from heeps.config import conf
from heeps.contrast.adi import adi
import multiprocessing as mpro
from functools import partial
from sys import platform
import time
import os

# email when finished
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_message'] = 'ADI simulation finished OK.'

# list of cases
cases = [{'band':'L', 'mode':'ELT', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'L', 'mode':'RAVC','pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'L', 'mode':'CVC', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'L', 'mode':'APP', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'L', 'mode':'CLC', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'M', 'mode':'ELT', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'M', 'mode':'RAVC','pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'M', 'mode':'CVC', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'M', 'mode':'APP', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'M', 'mode':'CLC', 'pscale': 5.21,'mag': 5,  'add_bckg':True},
         {'band':'N1','mode':'ELT', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'N1','mode':'CVC', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'N1','mode':'CLC', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'N2','mode':'ELT', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'N2','mode':'CVC', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'N2','mode':'CLC', 'pscale':10.78,'mag':-1.6,'add_bckg':True},
         {'band':'L', 'mode':'ELT', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'L', 'mode':'RAVC','pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'L', 'mode':'CVC', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'L', 'mode':'APP', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'L', 'mode':'CLC', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'M', 'mode':'ELT', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'M', 'mode':'RAVC','pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'M', 'mode':'CVC', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'M', 'mode':'APP', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'M', 'mode':'CLC', 'pscale': 5.21,'mag': 5,  'add_bckg':False},
         {'band':'N1','mode':'ELT', 'pscale':10.78,'mag':-1.6,'add_bckg':False},
         {'band':'N1','mode':'CVC', 'pscale':10.78,'mag':-1.6,'add_bckg':False},
         {'band':'N1','mode':'CLC', 'pscale':10.78,'mag':-1.6,'add_bckg':False},
         {'band':'N2','mode':'ELT', 'pscale':10.78,'mag':-1.6,'add_bckg':False},
         {'band':'N2','mode':'CVC', 'pscale':10.78,'mag':-1.6,'add_bckg':False},
         {'band':'N2','mode':'CLC', 'pscale':10.78,'mag':-1.6,'add_bckg':False}]

cases = [{'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':0},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':1},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':2},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':3},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':4},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':5},
         {'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':6}]
cases = [{'band':'L', 'mode':'RAVC','pscale':5.21, 'mag': 6,  'add_bckg':False, 'folder':10}]

# number of cases
ncases = len(cases)
# number of CPUs to use: one case per CPU
conf['cpucount'] = ncases

# define a wrapper function
def multi_cc(verbose, case):
    cube_duration = 3600
    cube_samp = 300
    adi_cube_duration = 3600
    adi_cube_samp = 300
    band = case['band']
    mode = case['mode']
    pscale = case['pscale']
    add_bckg = case['add_bckg']
    mag = case['mag']
    # optional sub-folder for onaxis PSFs
    folder = str(case.get('folder',''))
    tagname = '_%s'%folder
    path_offaxis = 'offaxis'
#    path_onaxis = 'cube_COMPASS_20181008_3600s_300ms_12000x256x256'
    path_onaxis = 'cube_COMPASS_20181008_3600s_300ms_12000x256x256_aberrations/%s/'%folder
    path_output = path_onaxis
    # call adi function
    adi(path_offaxis=path_offaxis, path_onaxis=path_onaxis, \
            path_output=path_output, verbose=verbose, \
            mode=mode, band=band, cube_duration=cube_duration, \
            cube_samp=cube_samp, adi_cube_duration=adi_cube_duration, \
            adi_cube_samp=adi_cube_samp, psc_simu=pscale, psc_inst=pscale, \
            add_bckg=add_bckg, mag=mag, tagname=tagname)
    # print stuff
    print('%s: Finished case %s.'%(time.strftime("%Y-%m-%d %H:%M:%S", \
            time.localtime()), folder))

# start ADI
print('%s: simulation starts, %s cases.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ncases))
if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
    print('Multiprocessing using %s cores'%conf['cpucount'])
    verbose = False
    func = partial(multi_cc, verbose)
    p = mpro.Pool(conf['cpucount'])
    p.map(func, cases)
    p.close()
    p.join()
else:
    for case in cases:
        verbose = True
        multi_cc(verbose, case)

# send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))
