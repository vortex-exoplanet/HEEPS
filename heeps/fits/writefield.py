from astropy.io import fits

def writefield(path, filename, field):

    fits.writeto(path + filename + '_r.fits', field.real, header=None, overwrite=True)
    fits.writeto(path + filename + '_i.fits', field.imag, header=None, overwrite=True)

    return



