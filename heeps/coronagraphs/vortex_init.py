import proper
import heeps.util.img_processing as impro
from heeps.coronagraphs.lens import lens
from ..fits import writefield
from ..fits import readfield
import numpy as np
import os.path

def vortex_init(conf):

    tmp_dir = conf['temp_dir']
    calib = conf['vortex_calib']
    my_file = str(tmp_dir+'zz_perf_'+calib+'_r.fits')
    
    if os.path.isfile(my_file):
        # read the pre-vortex field
        conf['psf_num'] = readfield(tmp_dir, 'zz_psf_%s'%calib)
        # read the theoretical vortex field
        conf['vvc'] = readfield(tmp_dir,'zz_vvc_%s'%calib) 
        # read the perfect-result vortex field
        conf['perf_num'] = readfield(tmp_dir,'zz_perf_%s'%calib)

    else:
        # load parameters
        lam = conf['lam']
        gridsize = conf['gridsize']
        beam_ratio = conf['beam_ratio']
        diam = conf['diam']
        f_lens = conf['focal']
        charge = conf['VC_charge']
        # create circular pupil
        wfo1 = proper.prop_begin(diam, lam, gridsize, beam_ratio)
        proper.prop_circular_aperture(wfo1, diam/2)
        proper.prop_define_entrance(wfo1)
        # propagate to vortex
        lens(wfo1, f_lens)
        # write the pre-vortex field
        conf['psf_num'] = wfo1.wfarr
        writefield(tmp_dir,'zz_psf_'+calib, conf['psf_num'])
        # vortex phase ramp is oversampled for a better discretization
        ramp_oversamp = 11.
        nramp = int(gridsize*ramp_oversamp)
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
        vvc = np.array(impro.resize_img(vvc_tmp.real, gridsize), dtype=complex)
        vvc.imag = impro.resize_img(vvc_tmp.imag, gridsize)
        phase_ramp = np.angle(vvc)
        # write the theoretical vortex field
        vvc_complex = np.array(np.zeros((gridsize, gridsize)), dtype=complex)
        vvc_complex.imag = phase_ramp
        conf['vvc'] = np.exp(vvc_complex)
        writefield(tmp_dir,'zz_vvc_'+calib, conf['vvc'])
        # apply vortex
        proper.prop_multiply(wfo1, conf['vvc'])
        # null the amplitude inside the Lyot Stop, and back propagate
        lens(wfo1, f_lens)
        proper.prop_circular_obscuration(wfo1, 1., NORM=True)
        lens(wfo1, -f_lens)
        # write the perfect-result vortex field
        conf['perf_num'] = wfo1.wfarr
        writefield(tmp_dir,'zz_perf_'+calib, conf['perf_num'])

    # shift the phase ramp
    conf['vvc'] = proper.prop_shift_center(conf['vvc'])

    return conf