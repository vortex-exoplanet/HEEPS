import heeps.util.img_processing as impro
from astropy.io import fits
import astropy.units as u
import numpy as np
import os
import time

# inputs
bands = ['L']#'L', 'M', 'N1', 'N2']
savename_fits = 'cc_raw_%s_%s_%s.fits'
# radius of CLC occulter in lam/D
r_CLC = 2.5
# modes per band
band_specs = {'L': {'modes': ['APP_adi', 'APP_raw'],#['RAVC', 'APP', 'CVC', 'CLC', 'ELT'],
                   'pscale': 5.47},
              'M': {'modes': ['RAVC', 'APP', 'CVC', 'CLC', 'ELT'],
                   'pscale': 5.47},
              'N1':{'modes': ['CVC', 'CLC', 'ELT'],
                   'pscale': 6.79},
              'N2':{'modes': ['CVC', 'CLC', 'ELT'],
                   'pscale': 6.79}}

# working repository
os.chdir(os.path.normpath(os.path.expandvars('$HOME/heeps_metis/output_files')))
cases = ['water_vapor/scao_only/yes',
           'water_vapor/scao_only/no',
           'water_vapor/scao_only/noTT']
cases = ['scao_only']

print('\n%s: producing raw contrast curves.'\
            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
for case in cases:
    print('%s'%case)

    for band in bands:
        print('%s'%band)
        pscale = band_specs[band]['pscale']
        modes = band_specs[band]['modes']
        for mode in modes:
            print('%s'%mode, end=', ')
            # off-axis PSF
            psf_OFF = fits.getdata(os.path.join(case, 'offaxis_PSF_%s_%s.fits'%(band, mode)))
            ndet = psf_OFF.shape[-1]
            # on-axis PSFs (cube)
            psf_ON = fits.getdata(os.path.join(case, 'onaxis_PSF_%s_%s.fits'%(band, mode)))
            psf_ON = np.array(psf_ON, ndmin=3)
            # average
            psf_ON_avg = np.mean(psf_ON, 0)
            # radial profiles
            rim = ndet // 2 # image radius
            y1 = impro.get_radial_profile(psf_OFF, (rim,rim), 1)[:-1]
            y2 = impro.get_radial_profile(psf_ON_avg, (rim,rim), 1)[:-1]
            # normalize by the peak of the off-axis PSF
            peak = np.max(y1)
            y1 /= peak
            y2 /= peak
            # x axis in lam/D
            x = pscale*1e-3*np.arange(rim)
            # CLC mode
            if 'CLC' in mode:
                mask = np.where(x>r_CLC)[0]
                x = x[mask]
                y1 = y1[mask]
                y2 = y2[mask]
            # save to hdu
            date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            header = fits.Header({'date':date, 'band':band, 'mode':mode})
            hdu = fits.PrimaryHDU(np.float32([x,y2]), header=header)
            hdu.writeto(os.path.join(case, savename_fits%(band, mode, case.replace('/','_'))), overwrite=True)
            print('saved to ', os.path.join(case, savename_fits%(band, mode, case.replace('/','_'))))