from .lens import lens
import heeps.util.img_processing as impro
from astropy.io import fits
import os.path
import numpy as np
import proper
from copy import deepcopy

def vortex_init(vortex_calib='', dir_temp='', diam_ext=37, lam=3.8, ngrid=1024, 
        beam_ratio=0.26, focal=660, vc_charge=2, verbose=False, **conf):

    '''
    
    Creates/writes vortex back-propagation fitsfiles, or loads them if files 
    already exist.
    The following parameters will be added to conf: 
        vortex_calib, psf_num, perf_num, vvc
    
    Returns: conf (updated and sorted)

    '''

    # update conf with local variables (remove unnecessary)
    conf.update(locals())
    [conf.pop(key) for key in ['conf', 'verbose'] if key in conf]

    # check if back-propagation params already loaded for this calib
    calib = '%s_%s_%3.4f'%(vc_charge, ngrid, beam_ratio)
    if vortex_calib == calib:
        return conf
        
    else:
        # check for existing file
        filename = os.path.join(dir_temp, 'calib_vvc_%s.fits'%calib)    
        if os.path.isfile(filename):
            if verbose is True:
                print('   loading vortex back-propagation params')
            data = fits.getdata(os.path.join(dir_temp, filename))
            # read the pre-vortex field
            psf_num = data[0] + 1j*data[1]
            # read the theoretical vortex field
            vvc = data[2] + 1j*data[3]
            # read the perfect-result vortex field
            perf_num = data[4] + 1j*data[5]

        # create files
        else:
            if verbose is True:
                print("   writing vortex back-propagation params")
            # create circular pupil
            wf_tmp = proper.prop_begin(diam_ext, lam, ngrid, beam_ratio)
            proper.prop_circular_aperture(wf_tmp, 1, NORM=True)
            # propagate to vortex
            lens(wf_tmp, focal)
            # pre-vortex field
            psf_num = deepcopy(wf_tmp.wfarr)
            # vortex phase ramp is oversampled for a better discretization
            ramp_oversamp = 11.
            nramp = int(ngrid*ramp_oversamp)
            start = -nramp/2 - int(ramp_oversamp)/2 + 0.5
            end   =  nramp/2 - int(ramp_oversamp)/2 + 0.5
            Vp = np.arange(start, end, 1.)
            # Pancharatnam Phase = arg<Vref,Vp> (horizontal input polarization)
            Vref = np.ones(Vp.shape)
            prod = np.outer(Vref, Vp)
            phiPan = np.angle(prod + 1j*prod.T)
            # vortex phase ramp exp(ilphi)
            ofst = 0
            ramp_sign = 1
            vvc_tmp = np.exp(1j*(ramp_sign*vc_charge*phiPan + ofst))
            vvc = np.array(impro.resize_img(vvc_tmp.real, ngrid), dtype=complex)
            vvc.imag = impro.resize_img(vvc_tmp.imag, ngrid)
            phase_ramp = np.angle(vvc)
            # theoretical vortex field
            vvc_complex = np.array(np.zeros((ngrid, ngrid)), dtype=complex)
            vvc_complex.imag = phase_ramp
            vvc = np.exp(vvc_complex)
            # apply vortex
            proper.prop_multiply(wf_tmp, vvc)
            # null the amplitude inside the Lyot Stop, and back propagate
            lens(wf_tmp, focal)
            proper.prop_circular_obscuration(wf_tmp, 1., NORM=True)
            lens(wf_tmp, -focal)
            # perfect-result vortex field
            perf_num = deepcopy(wf_tmp.wfarr)
            # write all fields 
            data = np.dstack((psf_num.real, psf_num.imag, vvc.real, vvc.imag,\
                perf_num.real, perf_num.imag)).T
            fits.writeto(os.path.join(dir_temp, filename), np.float32(data), overwrite=True)

        # shift the phase ramp
        vvc = proper.prop_shift_center(vvc)
        # add vortex back-propagation parameters at the end of conf
        conf = {k: v for k, v in sorted(conf.items())}
        conf.update(vortex_calib=calib, psf_num=psf_num, vvc=vvc, perf_num=perf_num)

        if verbose is True:
            print('   vc_charge=%s, ngrid=%s, beam_ratio=%3.4f'%\
                (vc_charge, ngrid, beam_ratio))

        return conf