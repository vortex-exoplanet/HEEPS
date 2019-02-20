#import matplotlib; matplotlib.use('agg') # to run on a headless server
from heeps.config import conf
from heeps.contrast.adi import adi
import multiprocessing as mpro
from functools import partial
from sys import platform
import time
import os

# email when finished
conf['send_to'] = None #'cdelacroix@uliege.be'
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

# number of CPUs to use: one case per CPU
cpucount = 3 #len(cases)+1

# define a wrapper function
def multi_cc(cases, ind):
    path_offaxis = '$HOME/INSTRUMENTS/METIS/offaxis'
    path_onaxis = '$HOME/INSTRUMENTS/METIS/cube_COMPASS_20181008_3600s_300ms_12000x256x256'
    cube_duration = 3600
    cube_samp = 300
    adi_cube_duration = 3600
    adi_cube_samp = 300
    band = cases[ind]['band']
    mode = cases[ind]['mode']
    pscale = cases[ind]['pscale']
    add_bckg = cases[ind]['add_bckg']
    mag = cases[ind]['mag']
    # call adi function
    adi(path_offaxis=path_offaxis, path_onaxis=path_onaxis, mode=mode, band=band,\
            cube_duration=cube_duration, cube_samp=cube_samp,\
            adi_cube_duration=adi_cube_duration, adi_cube_samp=adi_cube_samp,\
            psc_simu=pscale, psc_inst=pscale, add_bckg=add_bckg, mag=mag)
    # print stuff
    print('%s: Finished %s band, %s mode.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), band, mode))

# start ADI
print('%s: simulation starts, using %s cores.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cpucount))
if cpucount > 1:
    p = mpro.Pool(cpucount)
    func = partial(multi_cc, cases)
    p.map(func, range(len(cases)))
else:
    for ind in range(len(cases)):
        multi_cc(cases, ind)

# send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: " + '%s\n'%conf['send_message'], time.localtime()))
if conf['send_to'] is not None:
    os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
            conf['send_subject'], conf['send_to']))
