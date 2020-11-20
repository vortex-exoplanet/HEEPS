#!/usr/bin/env python3

import heeps
import numpy as np
from astropy.io import fits 
from multiprocessing import Pool
from functools import partial
import astropy.units as u

def get_disp_vals(ADCnu, steps):
    ''' Get the dispersion calculation '''
    filename = 'PSF_offsets_HCI-Llong_method_scale_zOpt{}.0.fits'.format(ADCnu)
    hdul = fits.open(conf['dir_input']+ filename)
    wavelength = hdul[0].data  # Sampled into 256 points
    zenith_angles = np.rad2deg(hdul[1].data) # Starts from 0 to 60 degrees (61 points)
    xoffsets = hdul[2].data
    yoffsets = hdul[3].data
    xdisp = hdul[4].data
    ydisp = hdul[5].data
    # Selecting HCI-L long band
    ind1, = np.where(wavelength==3.7)[0]
    ind2, = np.where(wavelength==3.9541178)[0]
#    L_hci_short = wavelength[132:170]  # 3.50 - 3.70µm
    L_hci_long = wavelength[ind1:ind2+1]  # 3.70 - 3.95µm
    ss = np.linspace(0, L_hci_long.shape[0]-1, steps, dtype=int)
    sel_L_long = L_hci_long[ss]
    fit_order = 5
    fit_res_x = np.zeros((steps, fit_order+1))
    fit_res_y = np.zeros((steps, fit_order+1))
    for i, index in enumerate(ss+ind1):
        fit_res_x[i] = np.polyfit(zenith_angles, xoffsets[index,:], fit_order)
        fit_res_y[i] = np.polyfit(zenith_angles, yoffsets[index,:], fit_order)

    fit_disp_y = np.zeros((steps, fit_order+1))
    for i, index in enumerate(ss+ind1):
        fit_disp_y[i] = np.polyfit(zenith_angles, ydisp[index,:], fit_order)        
    return sel_L_long, fit_res_x, fit_res_y, fit_disp_y

def get_transit_para():
    '''Get the transit parameters '''
    data = np.loadtxt(conf['dir_input']+ 'ha_parang_zangle_alfCen.txt')
    ha = data[:,0]
    para = data[:,1]
    star_zenith = data[:,2]
    return star_zenith

def disp_at_zenith_lambda(zenith, fitpara_x, fitpara_y, fit_disp_y, i, conf):
    if conf['adc'] == 'ADC':
        if conf['disp_tag'] == 'With_disp':
            tipx = np.polyval(fitpara_x[i], zenith)
            tipy = np.polyval(fitpara_y[i], zenith)
    if conf['adc'] == 'noADC':
        tipx = 0
        tipy = np.polyval(fit_disp_y[i], zenith)
    if conf['disp_tag'] == 'No_disp':
        tipx, tipy = 0, 0
    return tipx, tipy

conf = dict(
    dir_current = '/home/pacific/Desktop/eso_projects/metis_coronagraphs/compare_old_and_new_HEEPS_23_09_2019/HEEPS_20_11_2020/',
    bands = ['L'],
    modes = ['ELT'],
    cpu_count = None,
    add_scao = False
)

conf = heeps.config.read_config(**conf)

# sel_L_long = np.linspace(3.7, 3.95, 3)

# for conf['band'] in conf['bands']:
#     for conf['mode'] in conf['modes']:

def calc_single_poly_psf(sel_L_long, fitpara_x, fitpara_y, fit_disp_y, star_zeniths, conf, itr):
    zenith = star_zeniths[itr]
    phase_screens, pointing_errs, apo_drifts = heeps.wavefront.load_errors(verbose=True, **conf)
    steps = sel_L_long.shape[0]
    ppsf = np.zeros((steps, conf['ndet'], conf['ndet']), dtype=np.float32)
    for i in range(steps):
        lam = sel_L_long[i]*1e-6
        conf['band_specs']['L']['lam'] = lam
        conf = heeps.config.update_config(verbose=True, **conf) 

        wf = heeps.pupil.pupil(verbose=True, **conf)        

        dispx, dispy = disp_at_zenith_lambda(zenith, fitpara_x, fitpara_y, fit_disp_y, i, conf)
        t_res = (((lam/conf['diam_ext'])*u.rad).to(u.mas)).value
        # Calculates tip/tilt from mas to lambda/D
        dispx, dispy = dispx/t_res, dispy/t_res

        print('tip/tilt in lambda/D:', dispx, dispy)
        
        psf = heeps.wavefront.propagate_one(wf, phase_screen=phase_screens[itr], \
            tiptilt=[dispx, dispy], savefits=False, onaxis=True, **conf)
        ppsf[i] = psf
    return np.mean(ppsf, axis=0)

itr = 1
conf['disp_tag'] = 'With_disp'
conf['adc'] = 'ADC'
adc = 40 
steps = 5
conf['hfov'] = 1.5

conf = heeps.config.update_config(verbose=True, **conf) 

sel_L_long, fitpara_x, fitpara_y, fit_disp_y = get_disp_vals(adc, steps) 
star_zeniths = get_transit_para()[0:5]

# ppsf = calc_single_poly_psf(sel_L_long, fitpara_x, fitpara_y, fit_disp_y, star_zeniths, conf, itr)

p = Pool(4)
func = partial(calc_single_poly_psf, sel_L_long, fitpara_x, fitpara_y, fit_disp_y, star_zeniths, conf)
psf_zenith = np.array(p.map(func, range(star_zeniths.shape[0])), dtype=np.float32)
p.close()
p.join()

fits.writeto('test_no_scao.fits', psf_zenith, overwrite=True)

# star_zeniths = get_transit_para()

