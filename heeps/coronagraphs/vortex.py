import heeps.util.img_processing as impro
import numpy as np
import proper
import os
from ..fits import writefield 
from ..fits import readfield

def vortex(wfo, conf):
    
    # load parameters
    lam = conf['lam']
    gridsize = conf['gridsize']
    beam_ratio = conf['beam_ratio']
    diam = conf['diam']
    f_lens = conf['focal']
    charge = conf['VC_charge']
    tmp_dir = conf['temp_dir']
    
    # some more parameters
    
    if charge!=0:
        calib = str(charge)+str('_')+str(int(beam_ratio*100))+str('_')+str(gridsize)
        my_file = str(tmp_dir+'zz_perf_'+calib+'_r.fits')
        
        proper.prop_propagate(wfo, f_lens, 'inizio') # propagate wavefront
        proper.prop_lens(wfo, f_lens, 'focusing lens vortex') # propagate through a lens
        proper.prop_propagate(wfo, f_lens, 'CVC') # propagate wavefront
        
        if not os.path.isfile(my_file): # create the vortex for a perfectly circular pupil
            wfo1 = proper.prop_begin(diam, lam, gridsize, beam_ratio)
            proper.prop_circular_aperture(wfo1, diam/2)
            proper.prop_define_entrance(wfo1)
            proper.prop_propagate(wfo1, f_lens, 'inizio') # propagate wavefront
            proper.prop_lens(wfo1, f_lens, 'focusing lens vortex') # propagate through a lens
            proper.prop_propagate(wfo1, f_lens, 'CVC') # propagate wavefront     
            writefield(tmp_dir,'zz_psf_'+calib, wfo1.wfarr) # write the pre-vortex field
            
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
            # vortex phase ramp exp(ilphi) with optional offset
            ofst = 0
            vvc_tmp = np.exp(1j*(charge*phiPan + ofst))
            vvc = np.array(impro.resize_img(vvc_tmp.real, gridsize), dtype=complex)
            vvc.imag = impro.resize_img(vvc_tmp.imag, gridsize)
            phase_ramp = np.angle(vvc)
            
            vvc_complex = np.array(np.zeros((gridsize, gridsize)), dtype=complex)
            vvc_complex.imag = phase_ramp
            vvc = np.exp(vvc_complex)
            vvc_tmp = 0.
            writefield(tmp_dir,'zz_vvc_'+calib, vvc) # write the theoretical vortex field
            
            proper.prop_multiply(wfo1, vvc)
            proper.prop_propagate(wfo1, f_lens, 'OAP2')
            proper.prop_lens(wfo1, f_lens)
            proper.prop_propagate(wfo1, f_lens, 'forward to Lyot Stop')
            proper.prop_circular_obscuration(wfo1, 1., NORM=True) # null the amplitude iside the Lyot Stop
            proper.prop_propagate(wfo1, -f_lens) # back-propagation
            proper.prop_lens(wfo1, -f_lens)
            proper.prop_propagate(wfo1, -f_lens)
            writefield(tmp_dir,'zz_perf_'+calib, wfo1.wfarr) # write the perfect-result vortex field
        
        # read the theoretical vortex field
        vvc = readfield(tmp_dir,'zz_vvc_'+calib)
        vvc = proper.prop_shift_center(vvc)
        scale_psf = wfo._wfarr[0,0]
        # read the pre-vortex field
        psf_num = readfield(tmp_dir,'zz_psf_'+calib)
        psf0 = psf_num[0,0]
        psf_num = psf_num/psf0*scale_psf
        # read the perfect-result vortex field
        perf_num = readfield(tmp_dir,'zz_perf_'+calib) 
        perf_num = perf_num/psf0*scale_psf
        # the wavefront takes into account the real pupil with the perfect-result vortex field
        wfo._wfarr = (wfo._wfarr - psf_num)*vvc + perf_num 
        
        proper.prop_propagate(wfo, f_lens, "propagate to pupil reimaging lens")  
        proper.prop_lens(wfo, f_lens, "apply pupil reimaging lens")
        proper.prop_propagate(wfo, f_lens, "lyot stop")
            
    return wfo
