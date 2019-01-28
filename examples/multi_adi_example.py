#import matplotlib; matplotlib.use('agg') # to run on a headless server
from heeps.contrast.adi import adi
import multiprocessing as mpro
from functools import partial
from sys import platform
import time
import os

# email when finished
conf['send_to'] = 'cdelacroix@uliege.be'
conf['send_subject'] = 'fenrir'
conf['send_message'] = 'ADI simulation finished'

# list of cases
cases = [{'band': 'Lp', 'mode': 'ELT', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'Lp', 'mode': 'VC', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'Lp', 'mode': 'RAVC', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'Mp', 'mode': 'ELT', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'Mp', 'mode': 'VC', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'Mp', 'mode': 'RAVC', 'pscale': 5.21, 'add_bckg': True},
         {'band': 'N1', 'mode': 'ELT', 'pscale': 10.78, 'add_bckg': True},
         {'band': 'N1', 'mode': 'VC', 'pscale': 10.78, 'add_bckg': True},
         {'band': 'N2', 'mode': 'ELT', 'pscale': 10.78, 'add_bckg': True},
         {'band': 'N2', 'mode': 'VC', 'pscale': 10.78, 'add_bckg': True},
         {'band': 'Lp', 'mode': 'ELT', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'Lp', 'mode': 'VC', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'Lp', 'mode': 'RAVC', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'Mp', 'mode': 'ELT', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'Mp', 'mode': 'VC', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'Mp', 'mode': 'RAVC', 'pscale': 5.21, 'add_bckg': False},
         {'band': 'N1', 'mode': 'ELT', 'pscale': 10.78, 'add_bckg': False},
         {'band': 'N1', 'mode': 'VC', 'pscale': 10.78, 'add_bckg': False},
         {'band': 'N2', 'mode': 'ELT', 'pscale': 10.78, 'add_bckg': False},
         {'band': 'N2', 'mode': 'VC', 'pscale': 10.78, 'add_bckg': False}]

def multi_cc(cases, ind):
    path_offaxis = '/mnt/disk4tb/METIS/offaxis'
    path_onaxis = '/mnt/disk4tb/METIS/cube_COMPASS_20181008_3600s_100ms'
    cube_duration = 3600
    cube_samp = 100
    band = cases[ind]['band']
    mode = cases[ind]['mode']
    pscale = cases[ind]['pscale']
    add_bckg = cases[ind]['add_bckg']
    # call adi function
    adi(path_offaxis=path_offaxis,path_onaxis=path_onaxis,mode=mode,\
            cube_duration=cube_duration,cube_samp=cube_samp,band=band,\
            psc_simu=pscale, psc_inst=pscale, add_bckg=add_bckg)
    # print stuff
    print('%s: Finished %s band, %s mode.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), band, mode))

""" Multiprocessing """
# number of CPUs to use 
cpucount = len(cases) + 1
p = mpro.Pool(cpucount)
func = partial(multi_cc, cases)
print('%s: simulation starts, using %s cores.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cpucount))
p.map(func, range(len(cases)))

# Send email when simulation finished
print(time.strftime("%Y-%m-%d %H:%M:%S: Simulation finished OK.\n", time.localtime()))
os.system('echo "%s" | mail -s "%s" %s'%(conf['send_message'], \
        conf['send_subject'], conf['send_to']))
