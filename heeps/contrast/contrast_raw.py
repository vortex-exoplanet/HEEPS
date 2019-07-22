#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import os.path
from astropy.io import fits
import heeps.util.img_processing as impro
import numpy as np
import astropy.units as u
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import time

# inputs
bands = ['L', 'M', 'N1', 'N2']
diam = 37
r_CLC = 2.5 # radius of CLC occulter in lam/D
xbin = 1
rim = 400
(xo,yo) = (rim-1,rim-1)
linestyles = ['-','--']
thinline = 0.8
figsize = (12,4)#None#
title_pad = -20
markevery = 0.12
markersize = 4
xticks_log = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
atm = False
ylim = (1e-7, 1e0) if atm else (1e-8, 1e0)
APP_replaced = True
plot_offaxis = False
loglog = True
add_top_xticks = True
save2hdu = False
xlabel = r"Angular separation $[\lambda/D]$"
ylabel = r"Raw contrast"
title = r"%s band" +' (%s)'%('with SCAO residuals' if atm else 'perfect wavefront')
savename_png = 'cc_raw%s_%s'%('_atm' if atm else '', 'log' if loglog else 'lin') + '_%s.png'
savename_fits = 'cc_raw%s'%('_atm' if atm else '') + '_%s_%s.fits'


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
current_dir = '/Users/cdelacroix/INSTRUMENTS/METIS/heeps_analysis/'
path_output = os.path.join(current_dir, 'output_files')
path_offaxis = os.path.join(current_dir, 'offaxis')
if atm is True:
    path_onaxis = os.path.join(current_dir, 'cube_COMPASS_20181008_3600s_300ms_12000x512x512_averaged')
else:
    path_onaxis = os.path.join(current_dir, 'onaxis')

for band in bands:
    lam = band_specs[band]['lam']
    pscale = band_specs[band]['pscale']
    modes = band_specs[band]['modes']
    colors = band_specs[band]['colors']
    markers = band_specs[band]['markers']
    fig = plt.figure(figsize=figsize)
    for mode, color, marker in zip(modes, colors, markers):
        # if APP vertical band was replaced by horizontal one
        replaced = '_replaced' if mode == 'APP' and APP_replaced is True else ''
        
        # off-axis PSF
        psf_OFF = fits.getdata(os.path.join(path_offaxis, 'offaxis_PSF_%s_%s.fits' \
                %(band, mode)))
        # save PSF size
        npsf = psf_OFF.shape[1]
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
        # x axis in lam/D
        x = npsf/2/rim*pscale*1e-3/xbin*np.arange(rim)[:-1]/(lam/diam*u.rad).to('arcsec').value
        if mode == 'CLC':
            mask = np.where(x>r_CLC)[0]
            x = x[mask]
            y1 = y1[mask]
            y2 = y2[mask]
        # curves
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
        # save (x,y) for each mode in fits.hdu (optional)
        if save2hdu is True:
            date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            header = fits.Header({'xlabel':xlabel,'ylabel':ylabel,\
                    'date':date, 'wavelnth':lam, 'band':band, 'mode':mode})
            hdu = fits.PrimaryHDU((x,y2), header=header)
            hdu.writeto(savename_fits%(band, mode), overwrite=True)
            print('saved to ', savename_fits%(band, mode))
    # figure
    ax = fig.gca()
    if loglog is True:
        ax.set_xlim(left=1)
        ax.loglog()
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xticks([1,2,5,10,20,50])
    else:
        ax.set_xlim(left=0)
        ax.set_yscale('log')
    plt.grid()
    plt.grid(which='minor', linestyle=':')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title%band, pad=title_pad)
    plt.legend()
    plt.ylim(ylim)
    if add_top_xticks is True:
        ax1_xlim = np.array(ax.get_xlim())
        ax2_xlim = ax1_xlim*(lam/diam*u.rad).to('arcsec').value
        ax2 = ax.twiny()
        ax2.set_xlim(ax2_xlim)
        if loglog is True:
            ax2.set_xscale('log')
            xt2 = [x for x in xticks_log if ax2_xlim[0] <= x <= ax2_xlim[1]]
            ax2.set_xticks(xt2)
        str_ticks = [str(round(x,2))+'"' for x in ax2.get_xticks()]
        ax2.set_xticklabels(str_ticks)

    plt.show(block=False)
    plt.savefig(os.path.join(path_output, savename_png%band), dpi=300, transparent=True)
    plt.close()
