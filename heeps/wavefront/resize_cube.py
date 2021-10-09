from heeps.util.img_processing import resize_cube, pad_img
from astropy.io import fits
import numpy as np
import time
import os.path

#folder = '/mnt/disk4tb/METIS/METIS_COMPASS_CUBES/COMPASS_592x592'
folder = '/mnt/disk4tb/METIS/METIS_CBASIC_CUBES/'
cube = fits.getdata(os.path.join(folder, 'cube_COMPASS_20181008_3600s_300ms.fits'))#[:100]
mask = fits.getdata(os.path.join(folder, 'Telescope_Pupil.fits'))
npts = 592
new_size = 285
pad_resize = 39.9988/37
new_name = 'cube_COMPASS_20181008_3600s_300ms_0piston_meters_scao_only_%s.fits'%new_size

def is_odd(n):
    return bool(n % 2)

def oversamp(precision, start, end, stop=1e6):
    scale = end/start
    size = np.array([(x, x*scale) for x in np.arange(start, stop) \
            if x*scale % 1 < precision and is_odd(x) is is_odd(round(x*scale))])
    return size[0] if np.any(size) else 0

t0 = time.time()
cube[:,mask==0] = np.nan
cube = (cube.T - np.nanmean(cube,(1,2))).T # removing piston
cube *= 1e-6 # in meters


n1,n2 = oversamp(1e-2, npts, npts*39.9988/37)
print(npts, n1, n2, int(n2)/int(n1)*37)
cube = resize_cube(cube, int(n1), cpu_count=None, verbose=True)
cube = np.float32([pad_img(x, int(n2), np.nan) for x in cube])
cube = resize_cube(cube, new_size, cpu_count=None, verbose=True)

cube.shape

fits.writeto(os.path.join(folder, new_name), np.float32(cube), overwrite=True)
print(time.time() - t0)
