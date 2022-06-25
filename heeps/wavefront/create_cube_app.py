from heeps.util.multiCPU import multiCPU
import numpy as np
from astropy.io import fits


band = 'L'
ndet = 293                      # image size (fov 0.8")
cw = 39                         # PSF core width
modes = ['APPIMG', 'APPLMS']    # APP modes
stripe_width = [57, 29]         # bright stripe (IMG:57, LMS:29)
vertical_stripe = [True, False] # vertical stripe (or horizontal)
remove_stripe = True

''' 
reduce memory usage !!
pos -> raw contrast PSF
neg -> ADI contrast PSF
'''
def build_APP_old(ndet, cw, pos, neg):
        # merge pos and neg cubes
        psf_raw = pos
        cx = int(ndet/2)
        psf_raw[:cx,:] = neg[:cx,:]
        # apply scaling to the vertical middle band
        xmin = int((ndet-cw)/2)
        xmax = int((ndet+cw)/2)
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

def build_APP(ndet, cw, sw, vs, pos, neg):
        # merge pos and neg cubes
        psf_raw = pos
        psf_raw[neg<pos] = neg[neg<pos]
        # remove the bright stripe
        if remove_stripe is True:
            cmin = int((ndet-cw)/2)
            cmax = int((ndet+cw)/2)
            smin = int((ndet-sw)/2)
            smax = int((ndet+sw)/2)
            if vs is True:
                psf_raw[:cmin,smin:smax] = np.transpose(psf_raw[smin:smax,:cmin], (1, 0))
                psf_raw[cmax:,smin:smax] = np.transpose(psf_raw[smin:smax,cmax:], (1, 0))
            else:
                psf_raw[smin:smax,:cmin] = np.transpose(psf_raw[:cmin,smin:smax], (1, 0))
                psf_raw[smin:smax,cmax:] = np.transpose(psf_raw[cmax:,smin:smax], (1, 0))
        return psf_raw

for mode, sw, vs in zip(modes, stripe_width, vertical_stripe):
    # off-axis PSF can be either pos or neg
    filename = 'offaxis_PSF_%s_%s'%(band, mode) + '%s.fits'
    pos = fits.getdata(filename%'_pos')
    fits.writeto(filename%'', pos, overwrite=True)
    # build on-axis PSF (stitch pos and neg)
    filename = 'onaxis_PSF_%s_%s'%(band, mode) + '%s.fits'
    pos = np.array(fits.getdata(filename%'_pos'), ndmin=3)
    neg = np.array(fits.getdata(filename%'_neg'), ndmin=3)
    case = 'create %s band %s cube'%(band, mode)
    pos = multiCPU(build_APP, case=case,
                    posargs=[ndet, cw, sw, vs],
                    posvars=[pos, neg])
    print('saving ' + filename%'')
    fits.writeto(filename%'', pos, overwrite=True)
    # (pos, neg) = multiCPU(build_APP_old, case=case, multi_out=True,
    #                             posargs=[ndet, cw],
    #                             posvars=[pos, neg])
    # print('saving ' + filename%'_raw')
    # fits.writeto(filename%'raw', pos, overwrite=True)
    # print('saving ' + filename%'_adi')
    # fits.writeto(filename%'adi', neg, overwrite=True)
