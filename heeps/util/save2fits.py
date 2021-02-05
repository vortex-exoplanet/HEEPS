from astropy.io import fits 
import os
import numpy as np

def save2fits(data, name, dir_output='output_files', band='L', mode='RAVC', **conf):

    ''' Save data to a fits file. '''
    
    os.makedirs(dir_output, exist_ok=True)
    try:
        filename = '%s.fits'%name%(band, mode)
    except TypeError:
        filename = '%s_%s_%s.fits'%(name, band, mode)
    filename = os.path.join(dir_output, filename)
    fits.writeto(filename, np.float32(data), overwrite=True)

    return filename