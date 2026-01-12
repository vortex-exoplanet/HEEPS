from heeps.optics import vortex_init
from heeps.util.round2 import round2odd
from heeps.util.save2pkl import save2pkl
import astropy.units as u
import numpy as np
from pathlib import Path
from heeps.config.definition_lyotstops import COLD_STOPS

def update_config(band='L', band_specs={'L':{}}, mode='RAVC', lam=3.8e-6, 
        pupil_img_size=40, diam_ext=37, diam_int=11, ngrid=1024, pscale=5.47, 
        hfov=1, ravc_calc=False, ravc_t=0.8 ,ravc_r=0.6, saveconf=False, 
        verbose=False, **conf):
    """Updates some configuration parameters before running the end-to-end 
    simulation.
    
    The following parameters are updated to match the selected spectral band:
        lam, pscale, flux_star, flux_bckg, npupil, diam_norm, beam_ratio, ndet,
        hfov, hfov_lamD
    For the RAVC mode, the following parameters are also updated:
        r_obstr, ravc_t, ravc_r 
    
    Returns: conf (updated and sorted)

    """
    
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
    if 'VC' in mode:
        conf = vortex_init(verbose=verbose, **conf)

    if verbose is True:
        if mode in ['RAVC']:
            print('   ravc_calc=%s, ravc_t=%3.4f, ravc_r=%3.4f'
                %(ravc_calc, ravc_t, ravc_r))
        print('   npupil=%s, pscale=%.4f mas, lam=%3.4E m'%(npupil, pscale, lam))
        print('   hfov=%s arcsec (-> %s lam/D)'%(round(hfov, 2), round(hfov_lamD, 2)))
        print('   detector size (ndet)=%s (%s lam/D)'%(ndet, round(hfov_lamD*2, 2)))

    # Auto-selection of Lyot stop for METIS
    if conf['select_lyot'] == 'auto':
        print('\n   Auto-selection Lyot stop from definition')

        if band in ['L', 'M']:
            if conf['mode'] in ['RAVC']:
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) /\
                     COLD_STOPS['RLS-LM']['f_lyot_stop']

            elif conf['mode'] in ['CVC']:
                # see also backup stop CLS-LM-b
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) /\
                    COLD_STOPS['CLS-LM']['f_lyot_stop']
            elif conf['mode'] in ['CLC']:
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
                    COLD_STOPS['ULS-LM']['f_lyot_stop']
            elif conf['mode'] in ['ELT']:
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
                    COLD_STOPS['SPM-LM']['f_lyot_stop']
            else:
                print('   no auto-selected Lyot stop for mode=%s'%conf['mode'])
        elif band in ['N1', 'N2']:
            if conf['mode'] in ['CVC']:
                # see also backup stop CLS-N-b
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
                    COLD_STOPS['CLS-N']['f_lyot_stop']
            elif conf['mode'] in ['CLC']:
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
                    COLD_STOPS['ULS-N']['f_lyot_stop']
            elif conf['mode'] in ['ELT']:
                # see also backup stop SPM-N-b
                conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
                    COLD_STOPS['SPM-N']['f_lyot_stop']
        else:
                print('   no auto-selected Lyot stop for mode=%s'%conf['mode'])

        if not(Path(conf['f_lyot_stop']).is_file()):
            print('   [ WARNING ] no file for auto-select Lyot stop not found at %s'%conf['f_lyot_stop'])

    elif conf['select_lyot'] != '': # string is not empty
        print(f'\n   Selecting Lyot stop name {conf['select_lyot']}')
        if Path(conf['f_lyot_stop']).is_file():
            print(f' [ WARNING ] Lyot stop name is set via select_lyot and f_lyot_stop file ({conf['f_lyot_stop']})is also found. Using select_lyot.')

        conf['f_lyot_stop'] = Path(conf['f_lyot_stop']) / \
            COLD_STOPS[conf['select_lyot']]['f_lyot_stop']

    print('\n')
    return conf