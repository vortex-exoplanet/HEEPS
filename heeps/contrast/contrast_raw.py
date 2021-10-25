import heeps.util.img_processing as impro
from astropy.io import fits
import astropy.units as u
import numpy as np
import os
import time

# inputs
bands = ['L','N2']#'L', 'M', 'N1', 'N2']
xbin = 1
ndet = 403      # must be odd
rim = ndet // 2 # image radius
r_CLC = 2.5     # radius of CLC occulter in lam/D
APP_replaced = True
savename_fits = 'cc_raw_%s_%s_%s.fits'

# modes per band
band_specs = {'L': {'modes': ['RAVC'],# 'APP', 'CLC'],
                   'pscale': 5.47},
              'M': {'modes': ['RAVC'],#'ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                   'pscale': 5.47},
              'N1':{'modes': ['CVC'],#'ELT', 'CVC', 'CLC'],
                   'pscale': 6.79},
              'N2':{'modes': ['CVC'],
                   'pscale': 6.79}}

# working repository
os.chdir(os.path.normpath(os.path.expandvars('$HOME/heeps_metis/output_files/cfull')))#'cbasic'
folders = ['water_vapor/scao_only/G_0.4_nLF_10',
           'water_vapor/all_effects/G_0.4_nLF_10']
folders = ['water_vapor/scao_only/no',
           'water_vapor/scao_only/yes',
           'water_vapor/scao_only/noTT']
folders = ['all_effects_misseg_qacits', 'point_jitter', 'point_drift',
           'all_effects_misseg', 'apo_drift', 'talbot', 'all_ncpa']

print('\n%s: producing raw contrast curves.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
for folder in folders:
    print('%s'%folder)

    for band in bands:
        print('   %s'%band, end=', ')
        pscale = band_specs[band]['pscale']
        modes = band_specs[band]['modes']
        for mode in modes:
            print('      %s'%mode)
            # if APP vertical band was replaced by horizontal one
            replaced = '_replaced' if mode == 'APP' and APP_replaced is True else ''
            # off-axis PSF
            psf_OFF = fits.getdata(os.path.join(folder, 'offaxis_PSF_%s_%s.fits' \
                    %(band, mode)))
            # save PSF size
            npsf = psf_OFF.shape[1]
            # resample
            psf_OFF_ndet = impro.resize_img(psf_OFF, ndet)
            # on-axis PSFs (cube)
            psf_ON = fits.getdata(os.path.join(folder, 'onaxis_PSF_%s_%s%s.fits' \
                    %(band, mode, replaced)))
            if psf_ON.ndim != 3:
                psf_ON = np.array(psf_ON, ndmin=3)
            # average
            psf_ON_avg = np.mean(psf_ON, 0)
            # resample
            psf_ON_ndet = impro.resize_img(psf_ON_avg, ndet)
            # radial profiles
            y1 = impro.get_radial_profile(psf_OFF_ndet, (rim,rim), xbin)[:-1]
            y2 = impro.get_radial_profile(psf_ON_ndet, (rim,rim), xbin)[:-1]
            # normalize by the peak of the off-axis PSF
            peak = np.max(y1)
            y1 /= peak
            y2 /= peak

            # clip small/negative values to background level
            # bkgd = np.var(psf_ON_avg/peak)
            # y1[y1<ylim[0]] = np.min(y1[y1>bkgd])
            # y2[y2<ylim[0]] = np.min(y2[y2>bkgd])

            # x axis in lam/D
            x = npsf/ndet*pscale*1e-3*np.arange(rim)

            # if mode == 'CLC':
            #     mask = np.where(x>r_CLC)[0]
            #     x = x[mask]
            #     y1 = y1[mask]
            #     y2 = y2[mask]

            # save to hdu
            date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            header = fits.Header({'date':date, 'band':band, 'mode':mode})
            hdu = fits.PrimaryHDU(np.float32([x,y2]), header=header)
            hdu.writeto(os.path.join(folder, savename_fits%(band, mode, folder.replace('/','_'))), overwrite=True)
            print('saved to ', os.path.join(folder, savename_fits%(band, mode, folder.replace('/','_'))))
