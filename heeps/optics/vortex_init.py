from .lens import lens
import heeps.util.img_processing as impro
from astropy.io import fits
import os.path
import numpy as np
import proper
from copy import deepcopy

def readfield(path, filename):
    data_r = fits.getdata(os.path.join(path, '%s_r.fits'%filename))
    data_i = fits.getdata(os.path.join(path, '%s_i.fits'%filename))
    return(data_r + 1j*data_i)

def writefield(path, filename, field):
    fits.writeto(os.path.join(path, '%s_r.fits'%filename), field.real, header=None, overwrite=True)
    fits.writeto(os.path.join(path, '%s_i.fits'%filename), field.imag, header=None, overwrite=True)

def vortex_init(conf, calib, verbose=False):

    tmp_dir = conf['dir_temp']
    my_file = os.path.join(tmp_dir, 'zz_perf_%s_r.fits'%calib)
    
    if os.path.isfile(my_file):
        if verbose is True:
            print("   loading vortex back-propagation fitsfiles")
        # read the pre-vortex field
        conf['psf_num'] = readfield(tmp_dir, 'zz_psf_%s'%calib)
        # read the theoretical vortex field
        conf['vvc'] = readfield(tmp_dir,'zz_vvc_%s'%calib) 
        # read the perfect-result vortex field
        conf['perf_num'] = readfield(tmp_dir,'zz_perf_%s'%calib)

    else:
        if verbose is True:
            print("   creating/writing vortex back-propagation fitsfiles")
        # load parameters
        lam = conf['lam']
        ngrid = conf['ngrid']
        beam_ratio = conf['npupil']/conf['ngrid']*(conf['diam_ext']/conf['pupil_img_size'])
        diam = conf['diam_ext']
        f_lens = conf['focal']
        charge = conf['vc_charge']
        # create circular pupil
        wf_tmp = proper.prop_begin(diam, lam, ngrid, beam_ratio)
        proper.prop_circular_aperture(wf_tmp, 1, NORM=True)
#        proper.prop_define_entrance(wf_tmp)
        # propagate to vortex
        lens(wf_tmp, f_lens)
        # write the pre-vortex field
        conf['psf_num'] = deepcopy(wf_tmp.wfarr)
        writefield(tmp_dir,'zz_psf_%s'%calib, conf['psf_num'])
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
        vvc_tmp = np.exp(1j*(ramp_sign*charge*phiPan + ofst))
        vvc = np.array(impro.resize_img(vvc_tmp.real, ngrid), dtype=complex)
        vvc.imag = impro.resize_img(vvc_tmp.imag, ngrid)
        phase_ramp = np.angle(vvc)
        # write the theoretical vortex field
        vvc_complex = np.array(np.zeros((ngrid, ngrid)), dtype=complex)
        vvc_complex.imag = phase_ramp
        conf['vvc'] = np.exp(vvc_complex)
        writefield(tmp_dir,'zz_vvc_%s'%calib, conf['vvc'])
        # apply vortex
        proper.prop_multiply(wf_tmp, conf['vvc'])
        # null the amplitude inside the Lyot Stop, and back propagate
        lens(wf_tmp, f_lens)
        proper.prop_circular_obscuration(wf_tmp, 1., NORM=True)
        lens(wf_tmp, -f_lens)
        # write the perfect-result vortex field
        conf['perf_num'] = deepcopy(wf_tmp.wfarr)
        writefield(tmp_dir,'zz_perf_%s'%calib, conf['perf_num'])

    # shift the phase ramp
    conf['vvc'] = proper.prop_shift_center(conf['vvc'])

    return conf