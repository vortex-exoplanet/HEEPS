from heeps.coronagraphs.lens import lens
from heeps.coronagraphs.vortex_init import vortex_init

def vortex(wfo, conf):

    if conf['VC_charge'] != 0:

        # load vortex calibration files
        calib = '%s_%s_%s'%(conf['VC_charge'], int(conf['beam_ratio']*100), conf['gridsize'])
        if conf.get('vortex_calib') != calib:
            conf['vortex_calib'] = calib
            conf = vortex_init(conf)

        # propagate to vortex
        lens(wfo, conf['focal'])
        
        # apply vortex
        scale_psf = wfo._wfarr[0,0]/conf['psf_num'][0,0]
        wf_corr = (conf['psf_num']*conf['vvc'] - conf['perf_num'])*scale_psf
        wfo._wfarr = wfo._wfarr*conf['vvc'] - wf_corr

        import proper
        from astropy.io import fits
        fits.writeto('phi_after_vortex.fits', proper.prop_get_phase(wfo), overwrite=True)
        fits.writeto('amp_after_vortex.fits', proper.prop_get_amplitude(wfo), overwrite=True)
        
        # propagate to lyot stop
        lens(wfo, conf['focal'])
            
    return wfo