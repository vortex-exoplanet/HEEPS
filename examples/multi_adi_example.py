#import matplotlib; matplotlib.use('agg') # to run on a headless server
from heeps.adi import adi
from multiprocessing import Pool
import os.path

""" Multiprocessing """
# number of CPUs to use 
ncpu = 8
p = Pool(ncpu)

""" working repository """
#folder = '/mnt/disk4tb/METIS/heeps-analysis/'
folder = '$HOME/INSTRUMENTS/METIS/heeps-analysis/'
path_offaxis = 'offaxis'
path_onaxis = '../cube_COMPASS_20180223_600s_100ms'
path_output = 'output_files'
# absolute paths
folder = os.path.expandvars(folder)
path_offaxis = os.path.join(folder, path_offaxis)
path_onaxis = os.path.join(folder, path_onaxis)
path_output = os.path.join(folder, path_output)

""" run HEEPS """
adi.adi(path_offaxis=path_offaxis, path_onaxis=path_onaxis, path_output=path_output)

# TODO: add a loop on bands and modes