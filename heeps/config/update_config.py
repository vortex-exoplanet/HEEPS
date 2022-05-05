from heeps.optics import vortex_init
from heeps.util.round2 import round2odd
from heeps.util.save2pkl import save2pkl
import astropy.units as u
import numpy as np

def update_config(band='L', band_specs={'L':{}}, mode='RAVC', lam=3.8e-6, 
        pupil_img_size=40, diam_ext=37, diam_int=11, ngrid=1024, pscale=5.47, 
        hfov=1, ravc_calc=True, ravc_t=0.8 ,ravc_r=0.6, saveconf=False, 
        verbose=False, **conf):
    
    '''
    
    Update config parameters. The following parameters will be updated to match 
    the selected spectral band:
        lam, pscale, flux_star, flux_bckg, npupil, diam_norm, beam_ratio, ndet,
        hfov, hfov_lamD
    For the RAVC mode, the following parameters will be updated:
        r_obstr, ravc_t, ravc_r 
    
    Returns: conf (updated and sorted)

    '''
    
    if verbose is True:
        print('Simulation config: band=%s, mode=%s'%(band, mode))
        print('\u203e'*18)
    # update band specs
    if np.any(band_specs.get(band)):
        conf.update(band_specs.get(band))
    # calculate pupil size: must be odd for PROPER
    lam = conf.get('lam', lam)
    pscale = conf.get('pscale', pscale)
    lam_npupil = pupil_img_size*ngrid*pscale*u.mas.to('rad')
    npupil = round2odd(lam_npupil/lam)
    # recalculate wavelength based on npupil
    lam = lam_npupil/npupil
    del(lam_npupil)
    # normalized pupil image size (to pass to conf)
    diam_norm = pupil_img_size/diam_ext
    # calculate beam ratio (to pass to conf)
    beam_ratio = npupil/ngrid*diam_ext/pupil_img_size
    # final image size on detector
    ndet = round2odd(2*hfov/pscale*1e3)
    if not(ndet % 2) :
        ndet += 1
    hfov = ndet/2*pscale/1e3
    hfov_lamD = hfov*u.arcsec.to('rad')/(lam/diam_ext)
    # RAVC parameters (Mawet2013)
    if mode in ['RAVC'] and ravc_calc is True:
        r_obstr = diam_int/diam_ext
        ravc_t = 1 - (r_obstr**2 + r_obstr*np.sqrt(r_obstr**2 + 8))/4
        ravc_r = r_obstr/np.sqrt(1 - ravc_t) if diam_int > 0 else 0
    # update conf with local variables (remove unnecessary)
    conf.update(locals())
    [conf.pop(key) for key in ['conf', 'saveconf', 'verbose'] if key in conf]
    # sort alphabetically
    conf = {k: v for k, v in sorted(conf.items())}
    # save conf as pickle file
    if saveconf is True:
        save2pkl('conf', **conf)
    # load vortex back-propagation fitsfiles
    if mode in ['RAVC', 'CVC']:
        conf = vortex_init(verbose=verbose, **conf)

    if verbose is True:
        if mode in ['RAVC']:
            print('   ravc_calc=%s, ravc_t=%3.4f, ravc_r=%3.4f'
                %(ravc_calc, ravc_t, ravc_r))
        print('   npupil=%s, pscale=%.4f mas, lam=%3.4E m'%(npupil, pscale, lam))
        print('   hfov=%s arcsec (-> ndet=%s, %s lam/D)\n'
            %(round(hfov, 2), ndet, round(hfov_lamD, 2)))

    return conf