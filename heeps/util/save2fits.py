from astropy.io import fits
from astropy.wcs import WCS
import os
import numpy as np

def save2fits(data, name, dir_output='output_files', band='L', mode='RAVC', 
        pscale=5.47, dit=0.3, **conf):

    ''' Save data to a fits file. '''
    
    os.makedirs(dir_output, exist_ok=True)
    try:
        filename = '%s.fits'%name%(band, mode)
    except TypeError:
        filename = '%s_%s_%s.fits'%(name, band, mode)
    filename = os.path.join(dir_output, filename)
    fits.writeto(filename, np.float32(data), overwrite=True)

    # compatibility with Scopesim: provide a WCS for the PSF, with the pixel scale and DIT.
    if 'PSF' in name:
        hdr = fits.getheader(filename)
        naxis1, naxis2 = hdr['NAXIS1'], hdr['NAXIS2']
        wcs = WCS(naxis=3)
        wcs.wcs.ctype = ['LINEAR', 'LINEAR', 'LOCAL']
        wcs.wcs.crpix = [(naxis1 + 1) / 2, (naxis2 + 1) / 2, 1]
        wcs.wcs.crval = [0., 0., 0.]
        wcs.wcs.cdelt = [pscale*1e-3, pscale*1e-3, dit]
        wcs.wcs.cunit = ['arcsec', 'arcsec', 's']

        with fits.open(filename, 'update') as hdul:
            
            # add WCS to the FITS header
            hdul[0].header.update(wcs.to_header())

            # append all scalar parameters from the sorted conf dict
            conf.update(dir_output=dir_output, band=band, mode=mode, pscale=pscale, dit=dit)
            conf = {k: v for k, v in sorted(conf.items())}
            for key, value in conf.items():
                if type(value) not in (dict, list, np.ndarray):
                    hdul[0].header.append((key[:8], value, key))        

    return filename