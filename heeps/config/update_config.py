from heeps.optics import vortex_init
from heeps.util.round2odd import round2odd
import collections
import astropy.units as u
import numpy as np

def update_config(mode='RAVC', mode_specs={'RAVC':{}}, band='L', band_specs={'L':{}}, 
        lam=3.81e-6, pupil_img_size=39.9988, diam_ext=36.905, diam_int=11.213, 
        ngrid=1024, pscale=5.47, hfov=1, ravc_calc=True, ravc_t=0.76 ,ravc_r=0.62, 
        vc_charge=2, vortex_calib='', verbose=False, **conf):
    
    '''
    
    Update config parameters. ATTENTION, the following parameters 
    will be updated to match the selected spectral band and HCI mode:
        lam, pscale, flux_star, flux_bckg, ls_dRspi, ls_dRint, npupil, 
        ndet, ravc_t, ravc_r

    '''
    
    if verbose is True:
        print('Update config: mode=%s, band=%s'%(mode, band))
    # update mode and band specs
    conf.update(mode_specs.get(mode))
    conf.update(band_specs.get(band))
    # calculate pupil size: must be odd for PROPER
    lam = conf.get('lam', lam)
    pscale = conf.get('pscale', pscale)
    lam_npupil = pupil_img_size*ngrid*pscale*u.mas.to('rad')
    npupil = round2odd(lam_npupil/lam)
    # recalculate wavelength based on npupil
    lam = lam_npupil/npupil
    # final image size on detector
    ndet = round2odd(2*hfov/pscale*1e3)
    if not(ndet % 2) :
        ndet += 1
    hfov = ndet/2*pscale/1e3
    hfov_lamD = hfov*u.arcsec.to('rad')/(lam/diam_ext)
    # update conf
    conf.update(
        mode=mode,
        mode_specs=mode_specs,
        band=band,
        band_specs=band_specs,
        lam=lam,
        pupil_img_size=pupil_img_size,
        diam_ext=diam_ext,
        diam_int=diam_int,
        ngrid=ngrid,
        npupil=npupil,
        pscale=pscale,
        hfov=hfov,
        ndet=ndet,
        vc_charge=vc_charge
    )
    # vortex parameters
    if mode in ['RAVC', 'CVC']:
        # load vortex back-propagation fitsfiles
        beam_ratio = npupil/ngrid*(diam_ext/pupil_img_size)
        calib = '%s_%s_%3.4f'%(vc_charge, ngrid, beam_ratio)
        if vortex_calib != calib:
            vortex_calib = calib
            conf = vortex_init(conf, calib, verbose=verbose)
        # RAVC parameters (Mawet2013)
        if mode in ['RAVC'] and ravc_calc is True:
            r_obstr = diam_int/diam_ext
            ravc_t = 1 - (r_obstr**2 + r_obstr*np.sqrt(r_obstr**2 + 8))/4
            ravc_r = r_obstr/np.sqrt(1 - ravc_t)
    # update conf
    conf.update(
        ravc_calc=ravc_calc,
        ravc_t=ravc_t,
        ravc_r=ravc_r,
        vc_charge=vc_charge,
        vortex_calib=vortex_calib
    )

    if verbose is True:
        if mode in ['RAVC']:
            print('   ravc_calc=%s, ravc_t=%3.4f, ravc_r=%3.4f'\
                %(ravc_calc, ravc_t, ravc_r))
        print('   npupil=%s, pscale=%s mas, lam=%3.4E m'%(npupil, pscale, lam))
        print('   ndet=%s, hfov=%s arcsec (%s lam/D)'%(ndet, round(hfov,2), \
            round(hfov_lamD,2)))
        print('')

    # sort alphabetically
    conf = collections.OrderedDict(sorted(conf.items()))

    return conf
