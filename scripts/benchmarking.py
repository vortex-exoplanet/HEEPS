'''
Similar to the demo notebook.
Intended to use for profiling. See usage below.

USAGE:

    python -m cProfile -o profile.prof benchmarking.py

To inspect the generated profile:
    pyprof2calltree -k -i profile.prof -o callgraph.dot

or use

    python -m pstats profile.prof
    sort cumtime
    stats 30
'''
import cProfile
import pstats

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
print(os.environ["OPENBLAS_NUM_THREADS"])
print(os.environ.get("OPENBLAS_NUM_THREADS"))
# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["NUMEXPR_NUM_THREADS"] = "1"

import heeps

import threadpoolctl
threadpoolctl.threadpool_limits(limits=1, user_api='blas')
print(threadpoolctl.threadpool_info())
# from vip_hci.config.utils_conf import print_pool_map_counters


conf = dict(
    dir_current = '$HOME/heeps_metis',  # specify a directory
    dir_output='output_files_gox',
    f_phase = 'wavefront/COMPASS_201810_RandomWind_100screens_meters.fits',
    nframes = 100,                      # number of SCAO phase screens selected
    cpu_count = 10,                     # number of physical CPUs
    hfov=0.8
)

# conf = dict(
#     dir_current = '$HOME/heeps_metis',  # specify a directory
#     dir_output='output_files_gox/test_N2',
#     f_phase = 'wavefront/cfull/cube_Cfull_20220929_3600s_300ms_wv_ncpa_rms_600_lag_1_nzer_20_N2_119.fits',
#     nframes = 12000,                      # number of SCAO phase screens selected
#     nstep=10,
#     nframes_avg=10,
#     duration=3600,
#     dit=0.3,
#     cpu_count = 6,                     # number of physical CPUs
#     hfov=1.1
# )

conf = heeps.config.read_config(verbose=False, **conf)

conf = heeps.config.update_config(saveconf=True, verbose=True, **conf) 

wf = heeps.pupil.pupil(savefits=True, verbose=True, **conf)

# heeps.wavefront.propagate(wf, onaxis=False, avg=True, savefits=True, verbose=True, **conf)

# heeps.wavefront.propagate(wf, onaxis=True, savefits=True, verbose=True, **conf)

sep, raw = heeps.contrast.cc_raw(savefits=True, verbose=True, **conf)

# import yappi
# yappi.clear_stats()
# yappi.set_clock_type('cpu') # 'wall'
# yappi.start()


conf['cpu_count'] = 1   # best is to set cpu_count = 1 for cc_adi !!
sep1, adi1 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)

# yappi.stop()
# yappi.get_func_stats().print_all()
# # yappi.get_func_stats().print_all()
# yappi.get_thread_stats().print_all()

# def my_func():
#     sep1, adi1 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)

# cProfile.run('my_func()', 'copy_profile.prof')
# pstats.Stats('copy_profile.prof').sort_stats('cumtime').print_stats(2)

# # run your pipeline...
# print_pool_map_counters(limit=50)

# import matplotlib.pyplot as plt

# xlabel = 'Angular separation $[arcsec]$'
# ylabel_adi = '5-$\sigma$ sensitivity (contrast)'
# ylabel_raw = 'raw contrast'
# fig = plt.figure(figsize=(12, 6))
# fig.subplots(2, 1, sharex=True)
# fig.subplots_adjust(hspace=0.02)
# axes = fig.axes
# axes[0].set_ylim(1e-7, 1e-2)
# axes[1].set_ylim(1e-7, 9e-3)
# axes[1].set_xlabel(xlabel)
# for j, (ax, ylabel) in enumerate(zip(axes, [ylabel_raw, ylabel_adi])):
#     ax.set_ylabel(ylabel)
#     ax.grid(True), ax.grid(which='minor', linestyle=':')
#     ax.loglog()
#     ax.xaxis.set_major_formatter(plt.ScalarFormatter())
#     ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
#     ax.set_xlim(0.02, 0.75)
# axes[0].plot(sep, raw, 'C0', label='RAVC', marker='d', markevery=0.12, markersize=4)
# axes[0].legend()
# axes[0].set_title('HCI modes; scao only; star mag L = 6')
# axes[1].plot(sep1, adi1, 'C0', label='RAVC', marker='d', markevery=0.12, markersize=4)
# # axes[1].plot(sep2, adi2, ':C0', label='background', marker='d', markevery=0.12, markersize=4)
# axes[1].legend(ncol=2, loc='upper right')
