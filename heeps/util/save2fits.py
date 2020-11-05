from astropy.io import fits 
import os
import numpy as np

def save2fits(data, name, dir_output='output_files', band='L', mode='RAVC', **conf):

    ''' Save data to a fits file. '''
    
    os.makedirs(dir_output, exist_ok=True)
    filename = os.path.join(dir_output, '%s_%s_%s.fits'%(name, band, mode))
    fits.writeto(filename, np.float32(data), overwrite=True)

    return filename