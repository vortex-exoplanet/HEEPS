from astropy.io import fits 
import os.path
import numpy as np

def save2fits(data, name, dir_output='output_files', band='L', mode='RAVC', **conf):

    ''' Save data to a fits file. '''

    filename = os.path.join(dir_output, '%s_%s_%s.fits'%(name, band, mode))
    fits.writeto(filename, np.float32(data), overwrite=True)

    return filename