import numpy as np
from astropy.io import fits
import multiprocessing as mpro
from functools import partial
from sys import platform
import time


planes = ['PSF']        # select which plane to merge [PSF, LS]
#modes = ['APPv1','APPv2']
modes = ['APP']

#bands = ['L','M']       # select bandwidth
#deltas =  [21, 25]      # width of the bright central band
#ndets = [377, 473]               # image size: 512, 256

bands = ['L']
deltas = [19]
ndets = [403]


conf = {}
conf['cpucount'] = None
conf['nframes'] = 12000

def build_APP(plane, mode, band, dpx, npx, pos, neg):
        # merge pos and neg cubes
        psf_raw = pos
        cx = int(npx/2)
        psf_raw[:cx,:] = neg[:cx,:]

        # apply scaling to the vertical middle band
        suffix = ''
        xmin = int((npx-dpx)/2)
        xmax = int((npx+dpx)/2)
        dark_level = np.median(psf_raw)
        bright_level = np.median(np.hstack((psf_raw[xmin:xmax,:xmin],psf_raw[xmin:xmax,xmax:])))
        #bright_level = np.median(psf_raw[xmin:xmax,:])
        psf_adi = np.array(psf_raw)
        psf_adi[xmin:xmax,:xmin] *= (dark_level/bright_level)
        psf_adi[xmin:xmax,xmax:] *= (dark_level/bright_level)
        
        # replace the vertical band by the horizontal band
        psf_raw[xmin:xmax,:xmin] = np.transpose(psf_raw[:xmin,xmin:xmax], (1, 0))
        psf_raw[xmin:xmax,xmax:] = np.transpose(psf_raw[xmax:,xmin:xmax], (1, 0))

        return (psf_raw, psf_adi)

''' 
reduce memory usage !!
cube1 = pos -> raw contrast PSF
cube2 = neg -> ADI contrast PSF
'''
for onoff in ['onaxis']:
    for plane in planes:
        for band, delta, ndet in zip(bands, deltas, ndets):
            for mode in modes:
                # cube1 = postive PSF => will store raw contrast PSF
                cube1 = np.array(fits.getdata('%s_%s_%s_%s_pos.fits'%(onoff, plane,band,mode)), ndmin=3)
                # cube2 = negative PSF => will store ADI contrast PSF
                cube2 = np.array(fits.getdata('%s_%s_%s_%s_neg.fits'%(onoff, plane,band,mode)), ndmin=3)
                if conf['cpucount'] != 1 and platform in ['linux', 'linux2', 'darwin']:
                    if conf['cpucount'] == None:
                        conf['cpucount'] = mpro.cpu_count() - 1
                    print('%s: %s band, %s mode, using %s cores.'\
                            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                            band, mode, conf['cpucount']))
                    p = mpro.Pool(conf['cpucount'])
                    func = partial(build_APP, plane, mode, band, delta, ndet)
                    (cube1, cube2) = np.array(p.starmap(func, zip(cube1, cube2))).transpose(1, 0, 2, 3)
                    #cube1 = cubes[:,0,:,:]
                    #cube2 = cubes[:,1,:,:]
                    p.close()
                    p.join()

                # TODO: fix when cpu_count = 1
                else:
                    print('%s: %s band, %s mode, using 1 core.'\
                            %(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), \
                            band, mode))
                    cube1 = np.zeros((conf['nframes'], ndet, ndet))
                    cube2 = np.zeros((conf['nframes'], ndet, ndet))
                    for i, (pos, neg) in enumerate(zip(cube1, cube2)):
                        psf_raw, psf_adi = build_APP(plane, mode, band, delta, ndet, pos, neg)
                        cube1[i,:,:] = psf_raw      
                        cube2[i,:,:] = psf_adi

                filename = '%s_%s_%s_%s_%s.fits'%(onoff, plane, band, mode, 'raw')
                fits.writeto(filename, cube1, overwrite=True)
                print('saving %s'%filename)
                filename = '%s_%s_%s_%s_%s.fits'%(onoff, plane, band, mode, 'adi')
                fits.writeto(filename, cube2, overwrite=True)
                print('saving %s'%filename)
