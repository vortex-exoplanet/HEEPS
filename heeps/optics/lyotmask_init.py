from heeps.util.coord import polar_coord
from heeps.util.img_processing import resize_img, pad_img
import numpy as np
import os.path
from astropy.io import fits
import proper

def lyotmask_init(lyotmask_calib='', dir_temp='', clc_diam=80, pscale=5.47, 
        magnification=100, ngrid=1024, beam_ratio=0.26, verbose=False, **conf):

    '''
    
    Creates/writes classical lyot masks, or loads them if files 
    already exist.
    The following parameters will be added to conf: 
        lyotmask_calib, lyotmask
    
    Returns: conf (updated and sorted)

    '''

    # update conf with local variables (remove unnecessary)
    conf.update(locals())
    [conf.pop(key) for key in ['conf', 'verbose'] if key in conf]

    # check if lyot mask already loaded for this calib
    calib = 'lyotmask_%s_%s_%3.4f'%(clc_diam, ngrid, beam_ratio)
    if lyotmask_calib == calib:
        return conf
        
    else:
        # check for existing file
        filename = os.path.join(dir_temp, '%s.fits'%calib)    
        if os.path.isfile(filename):
            if verbose is True:
                print('   loading lyot mask')
            lyotmask = fits.getdata(os.path.join(dir_temp, filename))

        # create file
        else:
            if verbose is True:
                print("   writing lyot mask")
            # calculate the size (in pixels) of the lyot mask, 
            # rounded up to an odd value
            nmask = int(np.ceil(clc_diam/pscale))
            nmask += 1 - nmask % 2
            # create a magnified lyot mask
            rmask = clc_diam/pscale/nmask
            r, t = polar_coord(nmask*magnification)
            lyotmask = np.zeros(np.shape(r))
            lyotmask[r > rmask] = 1
            # resize and pad with ones to amtch ngrid
            lyotmask = pad_img(resize_img(lyotmask, nmask), ngrid, 1)
            # save as fits file
            fits.writeto(os.path.join(dir_temp, filename), np.float32(lyotmask), overwrite=True)

        # shift the lyotmask amplitude
        lyotmask = proper.prop_shift_center(lyotmask)
        # add lyotmask amplitude at the end of conf
        conf = {k: v for k, v in sorted(conf.items())}
        conf.update(lyotmask_calib=calib, lyotmask=lyotmask)

        if verbose is True:
            print('   clc_diam=%s, ngrid=%s, beam_ratio=%3.4f'%\
                (clc_diam, ngrid, beam_ratio))

        return conf