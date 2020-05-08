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
bands = ['N1']#'L', 'M', 'N1', 'N2']
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
save2hdu = True
xlabel = r"Angular separation $[\lambda/D]$"
ylabel = r"Raw contrast"
title = r"%s band" +' (%s)'%('with SCAO residuals' if atm else 'perfect wavefront')
savename_png = 'cc_raw%s_%s'%('_atm' if atm else '', 'log' if loglog else 'lin') + '_%s.png'
savename_fits = 'cc_raw%s'%('_atm' if atm else '') + '_%s_%s_%s.fits'


# modes per band
band_specs = {'L': {'lam': 3.8332E-06 ,
                    'mag': 5,
                 'pscale': 5.47,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'M': {'lam': 4.7773E-06 ,
                    'mag': 5,
                 'pscale': 5.47,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'N1':{'lam': 8.6606E-06 ,#8.7e-6,
                    'mag': -1.6,
                 'pscale': 6.79,#10.78,
                  'modes': ['CVC'],#'ELT', 'CVC', 'CLC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']},
              'N2':{'lam': 11.251E-06 ,
                    'mag': -1.6,
                 'pscale': 6.79,#10.78,
                  'modes': ['CVC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']}}

# working repository
print('\n%s: producing raw contrast curves.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ncases))
folders = ['scao_only','all_effects']
for folder in folders:
    path_onaxis = 'output_files/%s/'%folder
    path_offaxis = path_onaxis
    path_output = path_onaxis

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
                hdu.writeto(os.path.join(path_output, savename_fits%(band, mode, folder)), overwrite=True)
                print('saved to ', savename_fits%(band, mode, folder))
