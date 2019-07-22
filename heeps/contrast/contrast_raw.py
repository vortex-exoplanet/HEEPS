#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import os.path
from astropy.io import fits
import heeps.util.img_processing as impro
import numpy as np

# inputs
bands = ['L', 'M', 'N1', 'N2']
xbin = 1
rim = 400
(xo,yo) = (rim-1,rim-1)
ylim = (1e-7, 1e0)
linestyles = ['-','--']
thinline = 0.8
markevery = 30
markersize = 4
APP_replaced = True
plot_offaxis = False
figsize = (12,4)
suffix = '_atm_long'

# modes per band
band_specs = {'L': {'lam': 3.8e-6,
                    'mag': 5,
                 'pscale': 5.21,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'M': {'lam': 4.8e-6,
                    'mag': 5,
                 'pscale': 5.21,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'N1':{'lam': 8.7e-6,
                    'mag': -1.6,
                 'pscale': 10.78,
                  'modes': ['ELT', 'CVC', 'CLC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']},
              'N2':{'lam': 11.5e-6,
                    'mag': -1.6,
                 'pscale': 10.78,
                  'modes': ['ELT', 'CVC', 'CLC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']}}

# working repository
path_offaxis = '/Users/cdelacroix/INSTRUMENTS/METIS/heeps_analysis/offaxis'
#path_onaxis = '/Users/cdelacroix/INSTRUMENTS/METIS/heeps_analysis/onaxis'
path_onaxis = '/Users/cdelacroix/INSTRUMENTS/METIS/heeps_analysis/cube_COMPASS_20181008_3600s_300ms_12000x512x512_averaged'
path_output = '/Users/cdelacroix/INSTRUMENTS/METIS/heeps_analysis/output_files'


# saved file name
savename = 'cc_raw_%s' + '%s.png'%suffix


def savefig(fignum, band, ylim=ylim, savename=savename):
    plt.figure(fignum)
    plt.yscale('log')
    plt.grid()
    plt.grid(which='minor', linestyle=':')
    plt.xlabel("Angular separation [arcsec]")
    plt.ylabel(r"Raw contrast")
    plt.title(r"%s band"%band)
    plt.legend()
    plt.xlim(left=0)
    plt.ylim(ylim)
    plt.show(block=False)
    plt.savefig(os.path.join(path_output, savename%band), dpi=300, transparent=True)
    plt.close()

for band in bands:
    pscale = band_specs[band]['pscale']
    modes = band_specs[band]['modes']
    colors = band_specs[band]['colors']
    markers = band_specs[band]['markers']
    x = pscale*1e-3/xbin*np.arange(rim)[:-1]
    plt.figure(1, figsize=figsize)
    for mode, color, marker in zip(modes, colors, markers):
        # if APP vertical band was replaced by horizontal one
        replaced = '_replaced' if mode == 'APP' and APP_replaced is True else ''
        
        # off-axis PSF
        psf_OFF = fits.getdata(os.path.join(path_offaxis, 'offaxis_PSF_%s_%s.fits' \
                %(band, mode)))
        # resample
        psf_OFF_rim = impro.resize_img(psf_OFF, 2*rim)
        # on-axis PSFs (cube)
        psf_ON = fits.getdata(os.path.join(path_onaxis, 'onaxis_PSF_%s_%s%s.fits' \
                %(band, mode, replaced)))
        if psf_ON.ndim != 3:
            psf_ON = np.array(psf_ON, ndmin=3)
        # average
        psf_ON_avg = np.mean(psf_ON, 0)
        # resample
        psf_ON_rim = impro.resize_img(psf_ON_avg, 2*rim)
        # radial profiles
        y1 = impro.get_radial_profile(psf_OFF_rim, (xo,yo), xbin)[:-1]
        y2 = impro.get_radial_profile(psf_ON_rim, (xo,yo), xbin)[:-1]
        # normalize by the peak of the off-axis PSF
        peak = np.max(y1)
        y1 /= peak
        y2 /= peak
        # clip small/negative values to background level
        bkgd = np.var(psf_ON_avg/peak)
        y1[y1<ylim[0]] = np.min(y1[y1>bkgd])
        y2[y2<ylim[0]] = np.min(y2[y2>bkgd])
        
        # figures
        if plot_offaxis is True:
            plt.plot(x, y2, color=color, marker=marker, \
                    markersize=markersize, markevery=markevery, \
                    linestyle=linestyles[0], label='%s on-axis'%mode)
            plt.plot(x, y1, color=color, linewidth=thinline, \
                    marker=marker, markersize=markersize, markevery=markevery, \
                    linestyle=linestyles[1], label='%s off-axis'%mode)
        else:
            plt.plot(x, y2, color=color, marker=marker, \
                    markersize=markersize, markevery=markevery, \
                    linestyle=linestyles[0], label='%s'%mode)
    savefig(1, band)
