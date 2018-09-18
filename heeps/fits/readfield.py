import numpy as np
from astropy.io.fits import getdata


def readfield(path, filename):

    data_r, hdr = getdata(path + filename + '_r.fits', header=True)
    data_i = getdata(path + filename + '_i.fits')

    field = np.array(data_r, dtype=complex)
    field.imag = data_i


    return(field)
