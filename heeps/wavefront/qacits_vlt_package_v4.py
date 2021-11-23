# -*- coding: utf-8 -*-

"""
Contains the functions for computing tip-tilt estimates 
based on the QACITS method for Vortex coronagraph of charge 2.

Version 4: 2018-12-07
    - Documentation completed
    - new get_delta_i_exact function, based on photutils module routines
    - new QACITS estimator selection (between outer/inner/full)
"""

__author__ = 'E. Huby, Obs. de Paris'
__all__ = ['run_qacits_vlt', 'quadrant_tiptilt',
           'bin_images',
           'subimage',
           'get_delta_i', 'get_delta_i_exact',
           'circle_mask',
           'convert_to_cube',
           'get_psf_flux', 'get_psf_coordinates',
           'get_aperture_flux',
           'gauss_2D', 'fit_gauss_2D',
           'display_tiptilt_target', 'display_tiptilt_sequence']

import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import warp
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from scipy.stats import linregress
from configobj import ConfigObj
try:
    import photutils
    _exact_default_ = True
except:
    _exact_default_ = False

import warnings
warnings.simplefilter("error", OptimizeWarning)

def run_qacits_vlt(parameter_file, 
                   sci_cube, sci_dit, 
                   psf_cube, psf_dit, 
                   image_sampling,
                   vortex_center_yx, psf_center_yx, 
                   n_bin = 1,
                   model_calibration=False, calib_tt=None, 
                   force = None, exact=True,
                   disp_plots  = False, verbose=False) :
    """
    Main QACITS function for operation with ERIS vortex coronagraph.
    
    Parameters
    ----------
    parameter_file : string
        name and path to the parameter file containing the relevant parameter 
        values for QACITS.
    sci_cube : 2D or 3D array
        input science coronagraphic image cube used to estimate the tip-tilt 
        based on the QACITS method. 
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
    sci_dit : float
        Detector integration time of the science images.
    psf_cube : 3D or 3D array
        off-axis PSF image cube used to normalize the flux. 
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
        If a cube is given, the PSF flux is estimated by the mean flux over
        all PSF frames.
    psf_dit : float
        Detector integration time of the non corongraphic PSF image.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    vortex_center_yx : tuple of floats
        coordinates of the vortex center position in pixels.
        vortex_center_yx = (vortex_center_y, vortex_center_x).
    psf_center_yx : tuple of floats
        coordinates of the approximate PSF center position in pixels.
        psf_center_yx = (psf_center_y, psf_center_x).
    n_bin : integer
        number of tip-tilt estimates that will be returned.
        must be comprised between 0 and n_img.
        - if 0: no frame averaging, every estimates are returned in a 2D array
                of dimensions (n_img, 2).
        - if 1: average all frames before computing the tip-tilt estimate.
        - if n_bin=integer < n_img: will return n_bin estimates, i.e. frames are 
                averaged by bins of n_img/n_bin before computing the estimates.
    model_calibration : boolean
        If True, the QACITS model for ERIS (slopes and cubic coefficient) will
        be computed using the input sci_cube, psf_cube and calib_tt.
    calib_tt : 1D array
        when in model calibration mode, the true tip-tilt amplitude must be 
        provided in order to perform the fit of the models.
    force : string
        force the QACITS estimator to use only the inner (force='inner') or the
        outer (force='outer') estimator.
    exact : boolean
        if True, the photometry in the quadrant will be exact (requires the
        photutils module). Otherwise the photometry measurement is limited by
        the pixel sampling.
    disp_plots : boolean
        If True, display results in plots.
    verbose : boolean
        prints various things in the console, for debugging purposes.

    Returns
    -------
    tt_xy : 2D array
        Final tip-tilt estimates in a 2D array of dimensions (n_img, 2).
        tt_xy[:,0] and tt_xy[:,1] are the tip-tilt amplitudes along the x and y
        directions respectively.
    sci_cube_processed: 2D or 3D array
        processed image cube, including the binned frames.
    """
    
    ### HARD CODED PARAMETERS ###############################################
    psf_rad_lbdd = 1. # lambda/D
    
    ### PSF flux estimate ####################################################
    if len(psf_cube.shape) == 2:
        psf_cube = np.expand_dims(psf_cube, 0)
        n_psf = 1
    else :
        n_psf = psf_cube.shape[0]
    
    if len(sci_cube.shape) == 2:
        ny, nx = sci_cube.shape
        sci_cube = np.expand_dims(sci_cube, 0)
        n_sci = 1
    else :
        n_sci, ny, nx = sci_cube.shape
    
    psf_cube_flux = []
    for i,p in enumerate(psf_cube):
        disp_plots = (disp_plots and (i == (n_psf-1)))
        p_flux = get_psf_flux(p, psf_rad_lbdd, image_sampling,
                                cyx_guess=psf_center_yx,
                                display=disp_plots,
                                do_norm=False)
        psf_cube_flux.append(p_flux)
    psf_cube_flux = np.array(psf_cube_flux)
    mean_psf_flux = np.mean(psf_cube_flux)
    #-- scaling with integration times
    mean_psf_flux = mean_psf_flux * sci_dit / psf_dit

    ### Image binning ########################################################
    if ((n_bin == 1) and n_sci > 1):
        # all images averaged
        sci_cube_binned = np.expand_dims(np.mean(sci_cube, axis=0),0)
    elif ((n_bin == 0) and n_sci > 1):
        sci_cube_binned = sci_cube.copy()
    elif n_sci > 1:
        # bin the images
        sci_cube_binned = np.zeros((n_bin, ny, nx))
        bin_width = n_sci // n_bin
        i0 = n_sci - bin_width * n_bin
        for i in range(n_bin):
            sci_cube_binned[i]= np.mean(sci_cube[i0+bin_width*i:i0+bin_width*(i+1),:,:], axis=0)
    else :
        sci_cube_binned = sci_cube.copy()
        
    ### Tip-tilt estimate ####################################################
    #-- load QACITS params
    qacits_params = ConfigObj(parameter_file, unrepr=True)
    tiptilt_estimate = quadrant_tiptilt(qacits_params, 
                                        sci_cube_binned,vortex_center_yx, 
                                        mean_psf_flux, image_sampling,
                                        model_calibration=model_calibration,
                                        calib_tt=calib_tt, force=force,
                                        verbose=verbose, exact=exact)    
    
    ## DISPLAY ###################################################################
    if disp_plots is True:    
        #-- Conversion in polar coordinates
        tt_rt = np.zeros_like(tiptilt_estimate)
        tt_rt[:,0] = np.sqrt(tiptilt_estimate[:,0]**2+tiptilt_estimate[:,1]**2)
        tt_rt[:,1] = np.arctan2(tiptilt_estimate[:,1],tiptilt_estimate[:,0]) * 180./np.pi
        #-- some parameters for display
        subim_width_lbdd = 3.
        inner_rad_lbdd  = 1.9
        outer_rad_lbdd  = 2.3
        #-- Quadrant analysis results for the last image
        img_i = -1
        print('\n Quadrant Analysis results'+
              '\n x = {0:1.2f} l/D  ;  y = {1:1.2f} l/D'.format(tiptilt_estimate[img_i,0],tiptilt_estimate[img_i,1])+
              '\n R = {0:1.2f} l/D  ;  t = {1:1.2f} degr'.format(tt_rt[img_i,0],tt_rt[img_i,1]))
        fig1, ax1 = display_tiptilt_target(sci_cube_binned[img_i], image_sampling, vortex_center_yx,
                                           tiptilt_estimate[img_i], subim_width_lbdd=subim_width_lbdd, 
                                           fig_num=41, tt_lim=1., plot_title='QUADRANT ANALYSIS',
                                           img_circ_rad=(inner_rad_lbdd, outer_rad_lbdd), 
                                           tt_circ_rad=(.2,.4,.6,.8,1.))
        #-- plot all estimates
        fig2, ax2 = display_tiptilt_sequence(tiptilt_estimate, 
                                             tt_circ_rad=(.2,.4,.6,.8,1.), 
                                             fignum=42)
        plt.show(block=False)

    return tiptilt_estimate, sci_cube_binned

def quadrant_tiptilt_v8(qacits_params, img_cube, vortex_center_yx, psf_flux, 
        image_sampling, calib_tt=None, model_calibration=False, force=None, 
        exact=_exact_default_, verbose=False):
    """
    QACITS tip-tilt estimator method optimized for the VLT instruments (ERIS, NEAR).
    
     Parameters
    ----------
    qacits_params: dict
        dict object containing the relevant parameter values for QACITS.
    img_cube : 2D or 3D array
        input science coronagraphic image cube used to estimate the tip-tilt 
        based on the QACITS method. 
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
    vortex_center_yx : tuple of floats
        coordinates of the vortex center position in pixels.
        vortex_center_yx = (vortex_center_y, vortex_center_x).
    psf_flux : float
        Flux of the PSF integrated in an area of radius 1 lambda/D.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    calib_tt : 1D array
        when in model calibration mode, the true tip-tilt amplitude must be 
        provided in order to perform the fit of the models.
    model_calibration : boolean
        If True, runs the calibration mode using sci_cube, psf_cube and calib_tt.
    force : string
        force the QACITS estimator to use only the inner (force='inner') or the
        outer (force='outer') estimator.
    exact : boolean
        if True, the photometry in the quadrant will be exact (requires the
        photutils module). Otherwise the photometry measurement is limited by
        the pixel sampling.
    verbose : boolean
        prints various things in the console, for debugging purposes.
        
    Returns
    -------
    final_est_xy : ndarray
        final tip-tilt (x,y) estimate
    qacits_params : dict
        updated qacits parameters
    """
    ### Unpack the parameters (dependant on the instrument) #################
    radii = qacits_params['radii']
    ratio = qacits_params['ratio']

    ### Compute the differential intensities in the 3 regions #################
    if len(img_cube.shape) == 2:
        img_cube = np.expand_dims(img_cube, 0)
    n_img = img_cube.shape[0]
    img0 = img_cube[0]

    #-- Create the corresponding masks
    masks = {}
    for key in radii:
        masks[key]= (circle_mask(img0, radii[key][1] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0]) * 
                 (1.-circle_mask(img0, radii[key][0] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0])))
    
    #-- Compute the differential intensities
    all_dixy = {}
    for key in masks:
        if exact is True:
            all_dixy1 = get_delta_i_exact(img_cube/psf_flux, 
                    radii[key][1] * image_sampling,
                    cx=vortex_center_yx[1], cy=vortex_center_yx[0])
            if radii[key][0] !=0 :
                all_dixy0 = get_delta_i_exact(img_cube/psf_flux, 
                        radii[key][0] * image_sampling,
                        cx=vortex_center_yx[1], cy=vortex_center_yx[0])
            else :
                all_dixy0 = 0.
            all_dixy[key] = all_dixy1 - all_dixy0

    #return all_dixy

        else :
            all_dixy[key] = get_delta_i(img_cube*masks[key]/psf_flux, 
                    cx=vortex_center_yx[1], cy=vortex_center_yx[0])
    
    #-- Debias the diff. int. computed on full area from the linear component
    #   (estimated from the outer area)
    all_dixy['full'] += (ratio * all_dixy['outer'])

    # Transform to mod-arg (modulus-argument)
    all_di_mod = {}
    all_di_arg = {}
    for key in ['inner','outer','full']:
        all_di_mod[key] = np.sqrt(all_dixy[key][:,0]**2 + all_dixy[key][:,1]**2)
        all_di_arg[key] = np.arctan2(all_dixy[key][:,1], all_dixy[key][:,0])

    #***** model calibration mode ********************************************
    if model_calibration is True:
        tt_fit_lim = {'inner':(0., 0.1),'outer':(0., 0.5),'full':(0.2, 0.5)}
        colors = {'inner':[0.,0.3,.7],'outer':[.7,0.,0.3],'full':[0.,.7,0.5]}
        plt.figure(num=1, figsize=(12,9))
        plt.clf()
        fig, ax = plt.subplots(nrows=3,ncols=2,num=1)
        fig.subplots_adjust(hspace=0)
        coeffs = {}
        for i, region in enumerate(tt_fit_lim):
            ind_x = np.where((calib_tt>tt_fit_lim[region][0]) & 
                             (calib_tt<tt_fit_lim[region][1]))[0]
            x = calib_tt[ind_x]
            yy  = all_di_mod[region]
            y = yy[ind_x]
            if region == 'full':
                y  = np.abs(y)**(1/3) # full estimator 
            coeff = linregress(x,y).slope
            fit_coeff = calib_tt*coeff
            if region == 'full':
                coeff = coeff**3
                fit_coeff = fit_coeff**3
            ax[i,0].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,0].set_ylabel('Normalized Diff. Intensity')
            ax[i,0].plot(calib_tt, yy, 'o', color=colors[region], alpha=.9, markersize=2,
                       label=region+r' - r = {0:.1f} to {1:.1f} $\lambda/D$'
                       .format(radii[region][0], radii[region][1]))
            ax[i,0].plot(calib_tt, fit_coeff, color=colors[region], alpha=.6, linestyle='--', 
                       label=r'Fit coeff. [{0:.2f}-{1:.2f}] $\lambda/D$ = {2:.3f}'
                       .format(tt_fit_lim[region][0],tt_fit_lim[region][1],coeff))
            ax[i,0].grid(color='.8',linestyle='--')
            ax[i,0].set_xlim(0.,)
            ax[i,0].legend()
            error = ((yy - fit_coeff) / yy) * 100
            #error = ((yy - fit_coeff) / psf_flux) * 100
            ax[i,1].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,1].set_ylabel('Model Error [%]')
            ax[i,1].plot(calib_tt, error, 
                       'o', markersize=2, color=colors[region], alpha=.6)
            ax[i,1].set_ylim(-20., 20.)
            ax[i,1].set_xlim(0.,)
            ax[i,1].grid(color='.8',linestyle='--')
            coeffs[region] = coeff

        qacits_params['inner_slope'] = coeffs['inner']
        qacits_params['outer_slope'] = coeffs['outer']
        qacits_params['full_coeff'] = coeffs['full']
        print('\nModel calibration results:'+
              '\nInner slope = {0:.3f}\nOuter slope = {1:.3f}\nFull coeff  = {2:.3f}'
              .format(*coeffs.values()))
    
    ### Inverse the models ####################################################
    inner_slope = qacits_params['inner_slope']
    outer_slope = qacits_params['outer_slope']
    full_coeff = qacits_params['full_coeff']
    phase_tolerance = qacits_params['phase_tolerance']
    modul_tolerance = qacits_params['modul_tolerance']
    small_tt_regime = qacits_params['small_tt_regime']
    large_tt_regime = qacits_params['large_tt_regime']
    
    #-- inner region: linear
    inner_est      = np.zeros((n_img, 2))
    inner_est[:,0] = all_di_mod['inner'] / inner_slope
    inner_est[:,1] = all_di_arg['inner'] + np.pi
    #-- outer region: linear
    outer_est      = np.zeros((n_img, 2))
    outer_est[:,0] = all_di_mod['outer'] / outer_slope
    outer_est[:,1] = all_di_arg['outer']
    #-- full region: cubic
    full_est      = np.zeros((n_img, 2))
    full_est[:,0] = np.abs(all_di_mod['full']/full_coeff)**(1./3.)
    full_est[:,1] = all_di_arg['full']    

    #-- Estimator selection: 
    final_est = np.zeros((n_img, 2))
    test_output = np.zeros((n_img, 3))
    if force == 'inner':
        final_est = inner_est
    elif force == 'outer':
        final_est = outer_est
    elif force == 'full':
        final_est = full_est
    else :
        for i in range(n_img):
            # modulus to be trusted for choosing tt regime
            outer_modulus = outer_est[i,0]
            
            # build complex phasors
            inner_phasor = inner_est[i,0] * np.exp(1j * inner_est[i,1]) 
            outer_phasor = outer_est[i,0] * np.exp(1j * outer_est[i,1])
            full_phasor  = full_est[i,0]  * np.exp(1j * full_est[i,1])
            
            # test estimate agreement
            #-- phase agreement: IN/OUT
            test_phasor = np.exp(1j * inner_est[i,1]) * np.exp(-1j * outer_est[i,1])
            inout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            in_out_phase_agreement = (np.abs(inout_test_phase) < phase_tolerance/180*np.pi)
            #-- phase agreement: OUT/FULL
            test_phasor = np.exp(1j * full_est[i,1]) * np.exp(-1j * outer_est[i,1])
            fullout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            full_out_phase_agreement = (np.abs(fullout_test_phase) < phase_tolerance/180*np.pi)
            full_out_modul_agreement = (np.abs(full_est[i,0]-outer_est[i,0]) < full_est[i,0]*modul_tolerance)
            if verbose is True :
                print(i, inner_est[i,0], outer_est[i,0], full_est[i,0], in_out_phase_agreement, 
                        np.abs(test_phasor)*180./np.pi, in_out_phase_agreement)
                print('{0:03d} -- '.format(i) +
                       '\n \t IN   {0:.3f} l/D {1:.1f} deg'.format(inner_est[i,0],inner_est[i,1]*180/np.pi)+
                       '\n \t OUT  {0:.3f} l/D {1:.1f} deg'.format(outer_est[i,0],outer_est[i,1]*180/np.pi)+
                       '\n \t FULL {0:.3f} l/D {1:.1f} deg'.format(full_est[i,0],full_est[i,1]*180/np.pi))
            
            ### SMALL TIPTILT REGIME
            if outer_modulus < small_tt_regime :
                if verbose is True :
                    print('\n \t > IN-OUT phase agreement is {}'.format(in_out_phase_agreement))
                if in_out_phase_agreement == True:
                    meanphasor = (inner_phasor + outer_phasor)/2.
                    if verbose is True :
                        print('\t => small1: in+out')
                else :
                    meanphasor = outer_phasor 
                    if verbose is True :
                        print('\t => small2: out')

            ### LARGE TIPTILT REGIME
            else :
                if verbose is True :
                    print('\n \t > FULL-OUT phase agreement is {}'.format(full_out_phase_agreement)+
                    '\n \t > FULL-OUT modulus agreement is {}'.format(full_out_modul_agreement))
                if full_out_phase_agreement == True:
                    if full_out_modul_agreement == True:
                        meanphasor = ( outer_phasor + full_phasor )/2.
                        if verbose is True :
                            print('\t => large1: full+out')
                    else :
                        if full_est[i,0] < 1.:
                            meanphasor = full_phasor
                        else:
                            meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                        if verbose is True :
                            print('\t => large2: full')
                else :
                    if full_est[i,0] < 1.:
                        meanphasor = full_phasor
                    else:
                        meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                    if verbose is True :
                        print('\t => large3: full')

            if verbose is True:
                final_phase = np.arctan2(np.imag(meanphasor), np.real(meanphasor))*180./np.pi
                print('\t    final estimator mod = {0:.3f} l/D phase = {1:.1f} deg'.format(np.abs(meanphasor), final_phase))
            final_est[i,0] = np.abs(meanphasor)
            final_est[i,1] = np.arctan2(np.imag(meanphasor),np.real(meanphasor))
            test_output[i,:] = np.array([inout_test_phase,fullout_test_phase,np.abs(full_est[i,0]-outer_est[i,0])])
    
    ### Final estimator in X,Y ################################################
    final_est_xy = np.zeros_like(final_est)
    final_est_xy[:,0] = final_est[:,0] * np.cos(final_est[:,1])
    final_est_xy[:,1] = final_est[:,0] * np.sin(final_est[:,1])

    full_estimate_output = np.ndarray((n_img, 11))
    full_estimate_output[:,0:2] = final_est_xy
    full_estimate_output[:,2:4] = inner_est
    full_estimate_output[:,4:6] = outer_est
    full_estimate_output[:,6:8] = full_est
    full_estimate_output[:,8:] = test_output

    return full_estimate_output


def quadrant_tiptilt_v7b(qacits_params, img_cube, vortex_center_yx, 
                     psf_flux, image_sampling,
                     model_calibration = False, calib_tt = None, 
                     force=None, exact=_exact_default_, verbose=False):
    """
    QACITS tip-tilt estimator method optimized for the VLT instruments (ERIS, NEAR).
    
     Parameters
    ----------
    qacits_params: dict
        ConfigObj object containing the relevant parameter values for QACITS.
    img_cube : 2D or 3D array
        input science coronagraphic image cube used to estimate the tip-tilt 
        based on the QACITS method. 
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
    vortex_center_yx : tuple of floats
        coordinates of the vortex center position in pixels.
        vortex_center_yx = (vortex_center_y, vortex_center_x).
    psf_flux : float
        Flux of the PSF integrated in an area of radius 1 lambda/D.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    model_calibration : boolean
        If True, the QACITS model for ERIS (slopes and cubic coefficient) will
        be computed using the input sci_cube, psf_cube and calib_tt.
    calib_tt : 1D array
        when in model calibration mode, the true tip-tilt amplitude must be 
        provided in order to perform the fit of the models.
    force : string
        force the QACITS estimator to use only the inner (force='inner') or the
        outer (force='outer') estimator.
    exact : boolean
        if True, the photometry in the quadrant will be exact (requires the
        photutils module). Otherwise the photometry measurement is limited by
        the pixel sampling.
    verbose : boolean
        prints various things in the console, for debugging purposes.
        
    Returns
    -------
    full_estimate_output : ndarray
        the output array is of dimensions (n_img,11). For each image, there is 
        11 elements:
            - [,0:2] the first 2 corresponds to the final tip-tilt (x,y) estimate, 
            to be taken into account for correction.
            - [,2:4] 'inner estimate' based on the analysis of the inner part of the img,
            is given as (r, theta) format (amplitude and phase in radians)
            - [,4:6] 'outer estimate' based on the analysis of an annulus in the img,
            is given as (r, theta) format (amplitude and phase in radians)
            - [6:8] 'full estimate' based on the analysis of the full img,
            is given as (r, theta) format (amplitude and phase in radians)
            - [,8:11] 'test_diff' gives the phase and amplitude difference 
            between some estimators, solely for debug purposes
        Tip-tilt estimates are ndarray of dimensions (n_img, 2), with the second 
        axis corresponding to the x and y directions for the final estimate, 
        and to modulus and phase for the other estimates (r, theta).
    """
    ### Unpack the parameters (dependant on the instrument) ###
    radii = qacits_params['radii']
    inner_slope = qacits_params['inner_slope']
    outer_slope = qacits_params['outer_slope']
    full_coeff = qacits_params['full_coeff']
    ratio = qacits_params['ratio']
    phase_tolerance = qacits_params['phase_tolerance'] / 180. * np.pi
    modul_tolerance = qacits_params['modul_tolerance']
    small_tt_regime = qacits_params['small_tt_regime']
    large_tt_regime = qacits_params['large_tt_regime']

    #***** model calibration mode *********************************************
    if model_calibration is True:
        inner_slope = 1.
        outer_slope = 1.
        full_coeff  = 1.
        tt_fit_lim  = {'inner': (0., 0.1),
                       'outer': (0., 0.5),
                       'full':  (0.2, 0.5)}
    ###########################################################################
    
    if len(img_cube.shape) == 2:
        img_cube = np.expand_dims(img_cube, 0)
    n_img = img_cube.shape[0]
    
    ### Compute the differential intensities in the 3 regions #################
    img0 = img_cube[0]
    #-- Create the corresponding masks
    masks = {}
    for key in radii:
        masks[key]= (circle_mask(img0, radii[key][1] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0]) * 
                 (1.-circle_mask(img0, radii[key][0] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0])) )
    
    #-- Compute the differential intensities
    all_dix = {}
    all_diy = {}
    all_di_mod = {}
    all_di_arg = {}

    for key in masks:
        if exact is True:
            all_dixy1 = get_delta_i_exact(img_cube/psf_flux, radii[key][1] * image_sampling, 
                                   cx=vortex_center_yx[1], 
                                   cy=vortex_center_yx[0])
            if radii[key][0] !=0 :
                all_dixy0 = get_delta_i_exact(img_cube/psf_flux, radii[key][0] * image_sampling,
                                       cx=vortex_center_yx[1], 
                                       cy=vortex_center_yx[0])
            else :
                all_dixy0 = 0.
                
            all_dixy = all_dixy1 - all_dixy0
            
        else :
            all_dixy = get_delta_i(img_cube*masks[key]/psf_flux, 
                                   cx=vortex_center_yx[1], 
                                   cy=vortex_center_yx[0])
        all_dix[key] = all_dixy[:,0]
        all_diy[key] = all_dixy[:,1]
        
        all_di_mod[key] = np.sqrt(all_dix[key]**2+all_diy[key]**2)
        all_di_arg[key] = np.arctan2(all_diy[key], all_dix[key])

    #-- Debias the diff. int. computed on full area from the linear component
    #   (estimated from the outer area)
    all_dix['full'] += (ratio * all_dix['outer'])
    all_diy['full'] += (ratio * all_diy['outer'])
    all_di_mod['full'] = np.sqrt(all_dix['full']**2+all_diy['full']**2)
    all_di_arg['full'] = np.arctan2(all_diy['full'], all_dix['full'])
    
    ### Inverse the models ####################################################
    
    #-- inner region: linear
    inner_est      = np.zeros((n_img, 2))
    inner_est[:,0] = all_di_mod['inner'] / inner_slope
    inner_est[:,1] = all_di_arg['inner'] + np.pi
    #-- outer region: linear
    outer_est      = np.zeros((n_img, 2))
    outer_est[:,0] = all_di_mod['outer'] / outer_slope
    outer_est[:,1] = all_di_arg['outer']
    #-- full region: cubic
    full_est      = np.zeros((n_img, 2))
    full_est[:,0] = np.abs(all_di_mod['full']/full_coeff)**(1./3.)
    full_est[:,1] = all_di_arg['full']    

    #***** model calibration mode ********************************************
    if model_calibration is True:
        colors = {'inner':[0.,0.3,.7],'outer':[.7,0.,0.3],'full':[0.,.7,0.5]}
        slopes = {}
        model_order = {'inner':1,'outer':1,'full':3}
        plt.figure(num=1, figsize=(12,9))
        plt.clf()
        fig, ax = plt.subplots(nrows=3,ncols=2,num=1)
        fig.subplots_adjust(hspace=0) #wspace=0
        for i, region in enumerate(tt_fit_lim):
            ind_x = np.where((calib_tt>tt_fit_lim[region][0]) & 
                             (calib_tt<tt_fit_lim[region][1]))[0]
            x = calib_tt[ind_x]
            if   region == 'inner':
                yy  = inner_est[:,0]
            elif region == 'outer':
                yy  = outer_est[:,0]
            elif region == 'full':
                yy  = full_est[:,0]
            y = yy[ind_x]
            lin_params = linregress(x,y)
            slopes[region] = lin_params.slope
#            ax[i].set_title(region)
            ax[i,0].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,0].set_ylabel('Normalized Diff. Intensity')
            ax[i,0].plot(calib_tt, yy**model_order[region], 'o', 
                       color=colors[region], alpha=.9, markersize=2,
                       label=region+r' - r = {0:.1f} to {1:.1f} $\lambda/D$'
                       .format(radii[region][0], radii[region][1]))
            ax[i,0].plot(calib_tt, calib_tt**model_order[region]*slopes[region]**model_order[region], 
                       color=colors[region], alpha=.6, linestyle='--', 
                       label=r'Fit coeff. [{0:.2f}-{1:.2f}] $\lambda/D$ = {2:.3f}'
                       .format(tt_fit_lim[region][0],tt_fit_lim[region][1],slopes[region]**model_order[region]))
            ax[i,0].grid(color='.8',linestyle='--')
            ax[i,0].set_xlim(0.,)
            ax[i,0].legend()
            
            error = ((yy**model_order[region] - 
                      calib_tt**model_order[region]*slopes[region]**model_order[region]) /
                      yy**model_order[region]) * 100.
            ax[i,1].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,1].set_ylabel('Model Error [%]')
            ax[i,1].plot(calib_tt, error, 
                       'o', markersize=2, color=colors[region], alpha=.6)
            ax[i,1].set_ylim(-20., 20.)
            ax[i,1].set_xlim(0.,)
            ax[i,1].grid(color='.8',linestyle='--')
            
#            ax[i,2].imshow(img_cube[12]*(masks[region]), cmap='plasma', 
#                           vmin=np.min(img_cube[12]), vmax=np.max(img_cube[12]))
#            ax[i,2].set_yticklabels([])
#            ax[i,2].set_yticks([])
#            ax[i,2].set_xticklabels([])
#            ax[i,2].set_xticks([])
        
        inner_slope = slopes['inner']
        outer_slope = slopes['outer']
        full_coeff  = slopes['full']**3
        print('\nModel calibration results:'+
              '\nInner slope = {0:.3f}\nOuter slope = {1:.3f}\nFull coeff  = {2:.3f}'
              .format(inner_slope,outer_slope,full_coeff))
    # *************************************************************************
    
    #-- Estimator selection: 
    #file  = open(filelog, 'w')
    final_est = np.zeros((n_img, 2))
    test_output = np.zeros((n_img, 3))
    if force == 'inner':
        final_est = inner_est
    elif force == 'outer':
        final_est = outer_est
    elif force == 'full':
        final_est = full_est
    else :
        for i in range(n_img):
            #file.write('\n \t Frame {:3d}'.format(frame))
            
            # modulus to be trusted for choosing tt regime
            outer_modulus = outer_est[i,0] #outer_est[i,0]
            
            # build complex phasors
            inner_phasor = inner_est[i,0] * np.exp(1j * inner_est[i,1]) 
            outer_phasor = outer_est[i,0] * np.exp(1j * outer_est[i,1])
            full_phasor  = full_est[i,0]  * np.exp(1j * full_est[i,1])
            
            # test estimate agreement
            #-- phase agreement: IN/OUT
            test_phasor = np.exp(1j * inner_est[i,1]) * np.exp(-1j * outer_est[i,1])
            inout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            in_out_phase_agreement = (np.abs(inout_test_phase) < phase_tolerance)
            #-- phase agreement: OUT/FULL
            test_phasor = np.exp(1j * full_est[i,1]) * np.exp(-1j * outer_est[i,1])
            fullout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            
            full_out_phase_agreement = (np.abs(fullout_test_phase) < phase_tolerance)
            full_out_modul_agreement = (np.abs(full_est[i,0]-outer_est[i,0]) < full_est[i,0]*modul_tolerance)
            
            if verbose is True :
                print(i, inner_est[i,0], outer_est[i,0], full_est[i,0], in_out_phase_agreement, np.abs(test_phasor)*180./np.pi, in_out_phase_agreement)
                print('{0:03d} -- '.format(i) +
                       '\n \t IN   {0:.3f} l/D {1:.1f} deg'.format(inner_est[i,0],inner_est[i,1]*180./np.pi)+
                       '\n \t OUT  {0:.3f} l/D {1:.1f} deg'.format(outer_est[i,0],outer_est[i,1]*180./np.pi)+
                       '\n \t FULL {0:.3f} l/D {1:.1f} deg'.format(full_est[i,0],full_est[i,1]*180./np.pi))
                 #file.write('{0:03d} -- '.format(i) +
                 #      '\n \t IN   {0:.3f} l/D {1:.1f} deg'.format(inner_est[i,0],inner_est[i,1]*180./np.pi)+
                 #      '\n \t OUT  {0:.3f} l/D {1:.1f} deg'.format(outer_est[i,0],outer_est[i,1]*180./np.pi)+
                 #      '\n \t FULL {0:.3f} l/D {1:.1f} deg'.format(full_est[i,0],full_est[i,1]*180./np.pi))
            
            ### SMALL TIPTILT REGIME
            if outer_modulus < small_tt_regime :
                if verbose is True :
                    print('\n \t > IN-OUT phase agreement is {}'.format(in_out_phase_agreement))
                    #file.write('\n \t > IN-OUT phase agreement is {}'.format(in_out_phase_agreement))
                if in_out_phase_agreement == True:
                    meanphasor = (inner_phasor + outer_phasor)/2.
                    if verbose is True :
                        print('\t => small1: in+out')
                        #file.write('\t => small1: in+out')
                else :
                    meanphasor = outer_phasor 
                    if verbose is True :
                        print('\t => small2: out')
                        #file.write('\t => small2: out')
            ### LARGE TIPTILT REGIME
            else :
                if verbose is True :
                    print('\n \t > FULL-OUT phase agreement is {}'.format(full_out_phase_agreement)+
                    '\n \t > FULL-OUT modulus agreement is {}'.format(full_out_modul_agreement))
                    #file.write('\n \t > FULL-OUT phase agreement is {}'.format(full_out_phase_agreement)+
                    #'\n \t > FULL-OUT modulus agreement is {}'.format(full_out_modul_agreement))
                if full_out_phase_agreement == True:
                    if full_out_modul_agreement == True:
                        meanphasor = ( outer_phasor + full_phasor )/2.
                        if verbose is True :
                            print('\t => large1: full+out')
                            #file.write('\t => large1: full+out')
                    else :
                        if full_est[i,0] < 1.:
                            meanphasor = full_phasor
                        else:
                            meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                        if verbose is True :
                            print('\t => large2: full')
                            #file.write('\t => large2: full')
                else :
                    if full_est[i,0] < 1.:
                        meanphasor = full_phasor
                    else:
                        meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                    if verbose is True :
                        print('\t => large3: full')
                        #file.write('\t => large3: full')

            if verbose is True:
                final_phase = np.arctan2(np.imag(meanphasor), np.real(meanphasor))*180./np.pi
                print('\t    final estimator mod = {0:.3f} l/D phase = {1:.1f} deg'.format(np.abs(meanphasor), final_phase))
                #file.write('\n \t    final estimator mod = {0:.3f} l/D phase = {1:.1f} deg'.format(np.abs(meanphasor), final_phase))
            final_est[i,0] = np.abs(meanphasor)
            final_est[i,1] = np.arctan2(np.imag(meanphasor),np.real(meanphasor))
            test_output[i,:] = np.array([inout_test_phase,fullout_test_phase,np.abs(full_est[i,0]-outer_est[i,0])])
    
    #file.close()
    
    ### Final estimator in X,Y ################################################
    final_est_xy = np.zeros_like(final_est)
    final_est_xy[:,0] = final_est[:,0] * np.cos(final_est[:,1])
    final_est_xy[:,1] = final_est[:,0] * np.sin(final_est[:,1])
    
#    ### All estimators in R, Theta for DEBUG & LOG purposes ###################
#    all_est = {'inner_rt':  inner_est,
#               'outer_rt':  outer_est,
#               'full_rt':   full_est,
#               'test_diff': test_output}
#    
#    full_estimate_output = (final_est_xy, all_est)
    
    full_estimate_output = np.ndarray((n_img, 11))
    full_estimate_output[:,0:2] = final_est_xy
    full_estimate_output[:,2:4] = inner_est
    full_estimate_output[:,4:6] = outer_est
    full_estimate_output[:,6:8] = full_est
    full_estimate_output[:,8:] = test_output

    return full_estimate_output

def quadrant_tiptilt(qacits_params, img_cube, vortex_center_yx, 
                     psf_flux, image_sampling,
                     model_calibration = False, calib_tt = None, 
                     force=None, exact=True, verbose=False):
    """
    QACITS tip-tilt estimator method optimized for the ERIS instrument.
    
     Parameters
    ----------
    qacits_params: dict
        ConfigObj object containing the relevant parameter values for QACITS.
    img_cube : 2D or 3D array
        input science coronagraphic image cube used to estimate the tip-tilt 
        based on the QACITS method. 
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
    vortex_center_yx : tuple of floats
        coordinates of the vortex center position in pixels.
        vortex_center_yx = (vortex_center_y, vortex_center_x).
    psf_flux : float
        Flux of the PSF integrated in an area of radius 1 lambda/D.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    model_calibration : boolean
        If True, the QACITS model for ERIS (slopes and cubic coefficient) will
        be computed using the input sci_cube, psf_cube and calib_tt.
    calib_tt : 1D array
        when in model calibration mode, the true tip-tilt amplitude must be 
        provided in order to perform the fit of the models.
    force : string
        force the QACITS estimator to use only the inner (force='inner') or the
        outer (force='outer') estimator.
    exact : boolean
        if True, the photometry in the quadrant will be exact (requires the
        photutils module). Otherwise the photometry measurement is limited by
        the pixel sampling.
    verbose : boolean
        prints various things in the console, for debugging purposes.
        
    Returns
    -------
    final_est_xy : 2D array
        Final tip-tilt estimates in a 2D array of dimensions (n_img, 2). 
        final_est_xy[:,0] and final_est_xy[:,1] are the tip-tilt amplitudes 
        along the x and y directions respectively.
    """
    ### Unpack the parameters (dependant on the instrument) ###
    radii = qacits_params['radii']
    inner_slope = qacits_params['inner_slope']
    outer_slope = qacits_params['outer_slope']
    full_coeff = qacits_params['full_coeff']
    ratio = qacits_params['ratio']
    phase_tolerance = qacits_params['phase_tolerance'] * np.pi
    modul_tolerance = qacits_params['modul_tolerance']
    small_tt_regime = qacits_params['small_tt_regime']
    large_tt_regime = qacits_params['large_tt_regime']

    #***** model calibration mode *********************************************
    if model_calibration is True:
        inner_slope = 1.
        outer_slope = 1.
        full_coeff  = 1.
        tt_fit_lim  = {'inner': (0., 0.1),
                       'outer': (0., 0.5),
                       'full':  (0.2, 0.5)}
    ###########################################################################
    
    if len(img_cube.shape) == 2:
        img_cube = np.expand_dims(img_cube, 0)
    n_img = img_cube.shape[0]
    
    ### Compute the differential intensities in the 3 regions #################
    img0 = img_cube[0]
    #-- Create the corresponding masks
    masks = {}
    for key in radii:
        masks[key]= (circle_mask(img0, radii[key][1] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0]) * 
                 (1.-circle_mask(img0, radii[key][0] * image_sampling, 
                              cx=vortex_center_yx[1], cy=vortex_center_yx[0])) )
    
    #-- Compute the differential intensities
    all_dix = {}
    all_diy = {}
    all_di_mod = {}
    all_di_arg = {}

    for key in masks:
        if exact is True:
            all_dixy1 = get_delta_i_exact(img_cube/psf_flux, radii[key][1] * image_sampling, 
                                   cx=vortex_center_yx[1], 
                                   cy=vortex_center_yx[0])
            if radii[key][0] !=0 :
                all_dixy0 = get_delta_i_exact(img_cube/psf_flux, radii[key][0] * image_sampling,
                                       cx=vortex_center_yx[1], 
                                       cy=vortex_center_yx[0])
            else :
                all_dixy0 = 0.
                
            all_dixy = all_dixy1 - all_dixy0
            
        else :
            all_dixy = get_delta_i(img_cube*masks[key]/psf_flux, 
                                   cx=vortex_center_yx[1], 
                                   cy=vortex_center_yx[0])
        all_dix[key] = all_dixy[:,0]
        all_diy[key] = all_dixy[:,1]
        
        all_di_mod[key] = np.sqrt(all_dix[key]**2+all_diy[key]**2)
        all_di_arg[key] = np.arctan2(all_diy[key], all_dix[key])

    #-- Debias the diff. int. computed on full area from the linear component
    #   (estimated from the outer area)
    all_dix['full'] += (ratio * all_dix['outer'])
    all_diy['full'] += (ratio * all_diy['outer'])
    all_di_mod['full'] = np.sqrt(all_dix['full']**2+all_diy['full']**2)
    all_di_arg['full'] = np.arctan2(all_diy['full'], all_dix['full'])
    
    ### Inverse the models ####################################################
    
    #-- inner region: linear
    inner_est      = np.zeros((n_img, 2))
    inner_est[:,0] = all_di_mod['inner'] / inner_slope
    inner_est[:,1] = all_di_arg['inner'] + np.pi
    #-- outer region: linear
    outer_est      = np.zeros((n_img, 2))
    outer_est[:,0] = all_di_mod['outer'] / outer_slope
    outer_est[:,1] = all_di_arg['outer']
    #-- full region: cubic
    full_est      = np.zeros((n_img, 2))
    full_est[:,0] = np.abs(all_di_mod['full']/full_coeff)**(1./3.)
    full_est[:,1] = all_di_arg['full']    

    #***** model calibration mode ********************************************
    if model_calibration is True:
        colors = {'inner':[0.,0.3,.7],'outer':[.7,0.,0.3],'full':[0.,.7,0.5]}
        slopes = {}
        model_order = {'inner':1,'outer':1,'full':3}
        plt.figure(num=1, figsize=(12,9))
        plt.clf()
        fig, ax = plt.subplots(nrows=3,ncols=2,num=1)
        fig.subplots_adjust(hspace=0) #wspace=0
        for i, region in enumerate(tt_fit_lim):
            ind_x = np.where((calib_tt>tt_fit_lim[region][0]) & 
                             (calib_tt<tt_fit_lim[region][1]))[0]
            x = calib_tt[ind_x]
            if   region == 'inner':
                yy  = inner_est[:,0]
            elif region == 'outer':
                yy  = outer_est[:,0]
            elif region == 'full':
                yy  = full_est[:,0]
            y = yy[ind_x]
            lin_params = linregress(x,y)
            slopes[region] = lin_params.slope
#            ax[i].set_title(region)
            ax[i,0].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,0].set_ylabel('Normalized Diff. Intensity')
            ax[i,0].plot(calib_tt, yy**model_order[region], 'o', 
                       color=colors[region], alpha=.9, markersize=2,
                       label=region+r' - r = {0:.1f} to {1:.1f} $\lambda/D$'
                       .format(radii[region][0], radii[region][1]))
            ax[i,0].plot(calib_tt, calib_tt**model_order[region]*slopes[region]**model_order[region], 
                       color=colors[region], alpha=.6, linestyle='--', 
                       label=r'Fit coeff. [{0:.2f}-{1:.2f}] $\lambda/D$ = {2:.3f}'
                       .format(tt_fit_lim[region][0],tt_fit_lim[region][1],slopes[region]**model_order[region]))
            ax[i,0].grid(color='.8',linestyle='--')
            ax[i,0].set_xlim(0.,)
            ax[i,0].legend()
            
            error = ((yy**model_order[region] - 
                      calib_tt**model_order[region]*slopes[region]**model_order[region]) /
                      yy**model_order[region]) * 100.
            ax[i,1].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,1].set_ylabel('Model Error [%]')
            ax[i,1].plot(calib_tt, error, 
                       'o', markersize=2, color=colors[region], alpha=.6)
            ax[i,1].set_ylim(-20., 20.)
            ax[i,1].set_xlim(0.,)
            ax[i,1].grid(color='.8',linestyle='--')
            
#            ax[i,2].imshow(img_cube[12]*(masks[region]), cmap='plasma', 
#                           vmin=np.min(img_cube[12]), vmax=np.max(img_cube[12]))
#            ax[i,2].set_yticklabels([])
#            ax[i,2].set_yticks([])
#            ax[i,2].set_xticklabels([])
#            ax[i,2].set_xticks([])
        
        inner_slope = slopes['inner']
        outer_slope = slopes['outer']
        full_coeff  = slopes['full']**3
        print('\nModel calibration results:'+
              '\nInner slope = {0:.3f}\nOuter slope = {1:.3f}\nFull coeff  = {2:.3f}'
              .format(inner_slope,outer_slope,full_coeff))
    # *************************************************************************
    
    #-- Estimator selection: 
    final_est = np.zeros((n_img, 2))
    if force == 'inner':
        final_est = inner_est
    elif force == 'outer':
        final_est = outer_est
    else :
        for i in range(n_img):
            
            # modulus to be trusted for choosing tt regime
            modulus = outer_est[i,0]
            
            # build complex phasors
            phasor1 = inner_est[i,0] * np.exp(1j * inner_est[i,1]) 
            phasor2 = outer_est[i,0] * np.exp(1j * outer_est[i,1])
            phasor3 = full_est[i,0]  * np.exp(1j * full_est[i,1])
            
            # test estimate agreement
            test_phasor = np.exp(1j * inner_est[i,1]) * np.exp(-1j * outer_est[i,1])
            test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            in_out_phase_agreement = (np.abs(test_phase) < phase_tolerance)
            in_out_modul_agreement = (np.abs(full_est[i,0]-outer_est[i,0]) < modulus*modul_tolerance)
            
            if verbose is True :
                print(inner_est[i,0], outer_est[i,0], full_est[i,0], in_out_modul_agreement, np.abs(test_phase)*180./np.pi, in_out_phase_agreement)
           
            ## Select estimator according to phasor agreement and tt regime
            if in_out_phase_agreement :
                if modulus < small_tt_regime :    
                    meanphasor = ( phasor1 + phasor2 )/2.
                    if verbose is True :
                        print('small', np.abs(meanphasor))
                else : 
                    if in_out_modul_agreement == True:
                        meanphasor = ( phasor2 + phasor3 )/2.
                        if verbose is True :
                            print('intermediate', np.abs(meanphasor))
                    else :
                        meanphasor = 0.
                        if verbose is True :
                            print('null', np.abs(meanphasor))
            else:
                if modulus > large_tt_regime :
                    meanphasor = ( phasor2 + phasor3 )/2.
                    if verbose is True :
                        print('large', np.abs(meanphasor))
                elif in_out_modul_agreement == True : 
                    meanphasor = phasor3
                    if verbose is True :
                        print('large2', np.abs(meanphasor))
                else :
                    meanphasor = 0. 
                    if verbose is True :
                        print('null2', np.abs(meanphasor))
                        
            final_est[i,0] = np.abs(meanphasor)
            final_est[i,1] = np.arctan2(np.imag(meanphasor),np.real(meanphasor))
    
    ### Final estimator in X,Y ################################################
    final_est_xy = np.zeros_like(final_est)
    final_est_xy[:,0] = final_est[:,0] * np.cos(final_est[:,1])
    final_est_xy[:,1] = final_est[:,0] * np.sin(final_est[:,1])
    
    return final_est_xy


def bin_images(sci_cube, n_bin):
    """
    Returns a cube of images averaged by bins of n_bin images.
    
    Parameters
    ----------
    sci_cube : 2D or 3D array
        input science coronagraphic single image or image cube.
        Dimensions are (ny, nx) for 1 frame or (n_img, ny, nx) for a cube.
    n_bin : integer
        number of images in the returned cube.
        must be comprised between 0 and n_img.
        - if 0: no frame averaging, the returned cube is a copy of the input.
        - if 1: returns the average of the whole cube.
        - if n_bin=integer < n_img: the returned cube is made of n_bin images, 
            each being the average of n_img/n_bin images of the input sci_cube.
    
    Returns
    -------
    sci_cube_binned: 2D or 3D array
        binned image cube.
    """
    
    if len(sci_cube.shape) == 2:
        sci_cube = np.expand_dims(sci_cube, 0)
        n_sci = 1
        ny, nx = sci_cube.shape
    else :
        n_sci, ny, nx = sci_cube.shape
    
    if ((n_bin == 1) and n_sci > 1):
        # all images averaged
        sci_cube_binned = np.expand_dims(np.mean(sci_cube, axis=0),0)
    elif ((n_bin == 0) and n_sci > 1):
        sci_cube_binned = sci_cube.copy()
    elif n_sci > 1:
        # bin the images
        sci_cube_binned = np.zeros((n_bin, ny, nx))
        bin_width = n_sci // n_bin
        i0 = n_sci - bin_width * n_bin
        for i in range(n_bin):
            sci_cube_binned[i]= np.mean(sci_cube[i0+bin_width*i:i0+bin_width*(i+1),:,:], axis=0)
    else :
        sci_cube_binned = sci_cube.copy()
    
    return sci_cube_binned


def subimage(image, cx, cy, quadrant_width, full_output=False):
    """ 
    This routine creates a sub-image from a given image centered on the
    coordinates center given by (cx,cy) and of half width given by quadrant_width. 
    This sub-image is not necessarily a perfect square:
    the exact width in each direction depends on the parity of 2*cxy (cxy
    being cx or cy, and can be a fractional number). There are two different
    cases:
        - the center falls closer to a position between two pixels 
        (i.e. 2*cxy is odd), then the width is even (2*quadrant_width)
        - the center falls closer to a position in the middle of one pixel 
        (i.e. 2*cxy is even), then the width is odd (2*quadrant_width+1)

    Parameters
    ----------
    image : array_like
        input image for which the sub-image needs to be created.
    cx : float
        x position of the sub-image center [pix], can be fractional.
    cy : float
        y position of the sub-image center [pix], can be fractional.
    full_output : {True, False}, bool optional
        If True, the coordinates of the lower left pixels are returned.

    Returns
    -------
    subimage : array_like
        Subimage.
    """

    ny, nx = image.shape

    x1 = int(np.ceil(cx - quadrant_width))
    if x1 < 0 or x1 > nx :
        x1 = 0

    x2 = int(np.floor(cx + quadrant_width) + 1.)
    if x2 < 0 or x2 > nx :
        x2 = nx

    y1 = int(np.ceil(cy - quadrant_width))
    if y1 < 0 or y1 > ny :
        y1 = 0

    y2 = int(np.floor(cy + quadrant_width) + 1.)
    if y2 < 0 or y2 > ny :
        y2 = ny

    subimage = image[y1:y2, x1:x2]

    if full_output:
        return subimage, x1, y1
    else:
        return subimage


def get_delta_i(img_cube, cx=None, cy=None):
    """ 
    Computes the differential intensities along the x and y axes.

    Parameters
    ----------
    image : array_like
        Input 2D array.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.

    Returns
    -------
    delta_i : array_like
        2D element containing the differential intensities measured along the
        x and y axes.
    """
    
    if len(img_cube.shape) == 2:
        img_cube= np.expand_dims(img_cube,0)

    ny, nx = img_cube[0].shape

    if cx == None :
        cx = (nx-1.) / 2.
    if cy == None :
        cy = (ny-1) / 2.

    delta_i = []

    for img in img_cube:
        Sx = np.cumsum(np.sum(img, axis=0))
        Sy = np.cumsum(np.sum(img, axis=1))

        Ix = Sx[-1] - 2. * np.interp(cx, np.arange(nx)+.5, Sx)
        Iy = Sy[-1] - 2. * np.interp(cy, np.arange(ny)+.5, Sy)

        delta_i.append(np.array([Ix, Iy]))

    delta_i = np.array(delta_i)
    
    return delta_i

def get_delta_i_exact(img_cube, radius, cx=None, cy=None):
    """ 
    Computes the differential intensities along the x and y axes, based on the
    photmetry routines in the module photutils.

    Parameters
    ----------
    image : array_like
        Input 2D array.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.

    Returns
    -------
    delta_i : array_like
        2D element containing the differential intensities measured along the
        x and y axes.
    """
    
    if len(img_cube.shape) == 2:
        img_cube= np.expand_dims(img_cube,0)

    ny, nx = img_cube[0].shape

    if cx == None :
        cx = (nx-1.) / 2.
    if cy == None :
        cy = (ny-1) / 2.

    delta_i = []

    y, x = np.indices((ny,nx))
    cx1 = np.floor(cx)-1
    cy1 = np.floor(cy)-1
    
    aper = photutils.CircularAperture((cx, cy), radius)
    
    for img in img_cube:
         
        obj_flux = photutils.aperture_photometry(img*(x<=cx1), aper, method='exact')
        x_flux1 = obj_flux['aperture_sum'][0]
        obj_flux = photutils.aperture_photometry(img*(x<=(cx1+1)), aper, method='exact')
        x_flux2 = obj_flux['aperture_sum'][0]
        obj_flux = photutils.aperture_photometry(img, aper, method='exact')
        x_flux3 = obj_flux['aperture_sum'][0]
        Sx = [x_flux1,x_flux2,x_flux3]
        
        obj_flux = photutils.aperture_photometry(img*(y<=cy1), aper, method='exact')
        y_flux1 = obj_flux['aperture_sum'][0]
        obj_flux = photutils.aperture_photometry(img*(y<=(cy1+1)), aper, method='exact')
        y_flux2 = obj_flux['aperture_sum'][0]
        obj_flux = photutils.aperture_photometry(img, aper, method='exact')
        y_flux3 = obj_flux['aperture_sum'][0]
        Sy = [y_flux1,y_flux2,y_flux3]
                
        Ix = Sx[-1] - 2. * np.interp(cx, [cx1+.5,cx1+1.5,nx], Sx)
        Iy = Sy[-1] - 2. * np.interp(cy, [cy1+.5,cy1+1.5,ny], Sy)
        
        
        
        delta_i.append(np.array([Ix, Iy]))

    delta_i = np.array(delta_i)
    
    return delta_i


def circle_mask(image, radius_pix, cx=None, cy=None):
    """ 
    Creates a circular mask with values of 1 or 0, of the same dimensions as the
    input image, and defined by a radius and center coordinates. 
    Values are 1 inside the circle, 0 outside.

    Parameters
    ----------
    image : array_like
        Input 2D array.
    radius_pix : float
        Radius of the mask in pixels.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.

    Returns
    -------
    circle_mask : array_like
        2D array containing 1 and 0 values only.
    """
    ny, nx = image.shape

    if cx == None :
        cx = (nx-1.) / 2.
    if cy == None :
        cy = (ny-1) / 2.

    cx2 = np.round(cx*2.)/2.
    cy2 = np.round(cy*2.)/2.

    gridy, gridx = np.indices(image.shape)
    gridx = gridx - cx2
    gridy = gridy - cy2
    gridr = np.sqrt(gridx**2. + gridy**2.)

    circle_mask = gridr < radius_pix

    return circle_mask

def get_psf_flux_v7(psf_cube, radius_lbdd, image_sampling, cyx_guess,
                 t_int=1., display=False, do_norm = False, no_fit = False):
    """
    Estimates the flux in a circular aperture of the psf input image by 
    integrating the pixel counts in a circular area. 
    Integration time of the image can also be provided
    in order to normalize the final flux estimate.

    Parameters
    ----------
    psf_cube : array_like
        Off-axis PSF single image or cube of images.
        Dimensions should be (n_psf, nx, ny).
    radius_lbdd : float
        Radius of the mask in lambda/D.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    cyx_guess: tuple of floats
        approximate coordinates of the PSF center, used as a first guess of the
        2D Gaussian fit.
    t_int : float, optional
        Integration time of the psf image.
    do_norm : boolean
        if True, the flux is normalized such that the integrated flux of a flat
        image (1. everywhere) in an area of radius 1 lambda/D equals 1.
        Note: should be False for quadrant analysis
              should be True for annular analysis
    display : boolean
        If True, the off-axis PSF is displayed with a cross indicating the 
        fitted position of the 2D Gaussian profile.
    no_fit: boolean
        if True, the position of the PSF is not fitted, and the aperture 
        photometry is computed based on the input cyx_guess

    Returns
    -------
    psf_flux : float
        Flux of the psf integrated in the circular area and normalized by its
        integration time if provided.
    """
    
    if len(psf_cube.shape) == 2:
        psf_cube = np.expand_dims(psf_cube, 0)

    
    #### fit the position of the PSF
    # sub-image for fitting
    subim_width  = 4 * image_sampling # half sub-image width
    
    psf_cube_flux = [] 
    
    for i, psf in enumerate(psf_cube):
        psf_subimage = subimage(psf, cyx_guess[1], cyx_guess[0], 
                                subim_width, full_output=False)
        
        if no_fit is False :
            psf_cy, psf_cx = (get_psf_coordinates(psf_subimage, 
                                             cyx_guess = (subim_width, subim_width),
                                             sigma_guess = image_sampling/2.))[0]
        else :
            psf_cy, psf_cx = (subim_width, subim_width)
            
        radius_pix = radius_lbdd * image_sampling
        psf_flux   = get_aperture_flux(psf_subimage, radius_pix, 
                                       cx=psf_cx, cy=psf_cy, t_int=t_int)
        # Normalization (by the nb of pixels in an area of 1 lbd/D in radius)
        if do_norm is True:
            flat_flux  = get_aperture_flux(psf_subimage*0.+1., image_sampling, 
                                           cx=psf_cx, cy=psf_cy, t_int=t_int)
            psf_flux   = psf_flux / flat_flux

        if display is True:
            #        print('First guess 2D Gaussian params\n', first_guess)
            #        print('Fitted 2D Gaussian params\n', gauss_params)
            plt.figure(num=40, figsize=(4,4))
            plt.clf()
            plt.imshow(psf_subimage, cmap='plasma')
            plt.vlines(psf_cx, 0., 2*subim_width, color='w', linestyle='--')
            plt.hlines(psf_cy, 0., 2*subim_width, color='w', linestyle='--')
            thet = np.linspace(0., 2.*np.pi)
            plt.plot(psf_cx+radius_pix*np.cos(thet), 
                     psf_cy+radius_pix*np.sin(thet), 'w--')
            plt.title('2D Gaussian fit')
        
        psf_cube_flux.append(psf_flux)
    psf_cube_flux = np.array(psf_cube_flux)
    
    return psf_cube_flux

def get_psf_flux(psf_cube, radius_lbdd, image_sampling, cyx_guess,
                 t_int=1., display=False, do_norm = False):
    """
    Estimates the flux in a circular aperture of the psf input image by 
    integrating the pixel counts in a circular area. 
    Integration time of the image can also be provided
    in order to normalize the final flux estimate.

    Parameters
    ----------
    psf_cube : array_like
        Off-axis PSF single image or cube of images.
        Dimensions should be (n_psf, nx, ny).
    radius_lbdd : float
        Radius of the mask in lambda/D.
    image_sampling : float
        Image sampling given in pixels per lambda/D.
    cyx_guess: tuple of floats
        approximate coordinates of the PSF center, used as a first guess of the
        2D Gaussian fit.
    t_int : float, optional
        Integration time of the psf image.
    do_norm : boolean
        if True, the flux is normalized such that the integrated flux of a flat
        image (1. everywhere) in an area of radius 1 lambda/D equals 1.
        Note: should be False for quadrant analysis
              should be True for annular analysis
    display : boolean
        If True, the off-axis PSF is displayed with a cross indicating the 
        fitted position of the 2D Gaussian profile.

    Returns
    -------
    psf_flux : float
        Flux of the psf integrated in the circular area and normalized by its
        integration time if provided.
    """
    
    if len(psf_cube.shape) == 2:
        psf_cube = np.expand_dims(psf_cube, 0)

    
    #### fit the position of the PSF
    # sub-image for fitting
    subim_width  = 4 * image_sampling # half sub-image width
    
    psf_cube_flux = [] 
    
    for i, psf in enumerate(psf_cube):
        psf_subimage = subimage(psf, cyx_guess[1], cyx_guess[0], 
                                subim_width, full_output=False)
        
        psf_cy, psf_cx = (get_psf_coordinates(psf_subimage, 
                                             cyx_guess = (subim_width, subim_width),
                                             sigma_guess = image_sampling/2.))[0]
    
        radius_pix = radius_lbdd * image_sampling
        psf_flux   = get_aperture_flux(psf_subimage, radius_pix, 
                                       cx=psf_cx, cy=psf_cy, t_int=t_int)
        # Normalization (by the nb of pixels in an area of 1 lbd/D in radius)
        if do_norm is True:
            flat_flux  = get_aperture_flux(psf_subimage*0.+1., image_sampling, 
                                           cx=psf_cx, cy=psf_cy, t_int=t_int)
            psf_flux   = psf_flux / flat_flux

        if display is True:
            #        print('First guess 2D Gaussian params\n', first_guess)
            #        print('Fitted 2D Gaussian params\n', gauss_params)
            plt.figure(num=40, figsize=(4,4))
            plt.clf()
            plt.imshow(psf_subimage, cmap='plasma')
            plt.vlines(psf_cx, 0., 2*subim_width, color='w', linestyle='--')
            plt.hlines(psf_cy, 0., 2*subim_width, color='w', linestyle='--')
            thet = np.linspace(0., 2.*np.pi)
            plt.plot(psf_cx+radius_pix*np.cos(thet), 
                     psf_cy+radius_pix*np.sin(thet), 'w--')
            plt.title('2D Gaussian fit')
        
        psf_cube_flux.append(psf_flux)
    psf_cube_flux = np.array(psf_cube_flux)
    
    return psf_cube_flux

def convert_to_cube(img):
    """
    Checks if the input img is a cube or a single frame. If in the later case, 
    adds an extra dimension to to make it the same format as image cubes. 
    If a cube is provided, the routine has no effect.
    
    Parameters
    ----------
    img: 2D or 3D array
        input single image or image cube.
        
    Returns
    -------
    img_cube: 3D array
        image cube.
    """
    
    if len(img.shape) == 2:
        img_cube = np.expand_dims(img, 0)
        return img_cube
    else:
        return img

def get_psf_coordinates(psf_cube, cyx_guess=None, sigma_guess = None):
    """
    Estimates the coordinates of the PSF center in psf_cube.
    
    Parameters
    ----------
    psf_cube : array_like
        Off-axis PSF single image or cube of images.
        Dimensions should be (n_psf, nx, ny).
    cyx_guess: tuple of floats
        approximate coordinates of the PSF center, used as a first guess of the
        2D Gaussian fit.
    sigma_guess: float
        approximate sigmat parameter for the 2D Gaussian profile to be fitted.
    """
    
    psf_cube = convert_to_cube(psf_cube)
    n_img, ny, nx = psf_cube.shape
    psf_coord = np.zeros((n_img,2))
    
    if cyx_guess is None:
        cyx_guess = (ny/2.,  nx/2.)
    
    if sigma_guess is None:
        sigma_guess = 5.
    
    for i, psf in enumerate(psf_cube):
        # first guess parameters
        init_yx_max = np.unravel_index(psf.argmax(), 
                                       psf.shape)
        first_guess_1 = [psf.max(),                # amplitude
                         init_yx_max[1], init_yx_max[0],      # cx, cy
                         sigma_guess, sigma_guess,# sigma_x, sigma_y
                         0., 0.]                              # theta, offset
        
        first_guess_2 = [psf.max(),                # amplitude
                         cyx_guess[0], cyx_guess[1],            #cy, cx
                         sigma_guess, sigma_guess,# sigma_x, sigma_y
                         0., 0.]                              # theta, offset
        
        # fit the 2D Gaussian
        try :
            first_guess = first_guess_2
            #print('psf shape',psf.shape)
            #print('first guess',first_guess)
            gauss_params = fit_gauss_2D(psf, first_guess = first_guess)
            # fitted center of the Gaussian
            psf_cx = gauss_params[1]
            psf_cy = gauss_params[2]
        except RuntimeError:
            try :
                print('WARNING: first Gaussian fit failed.')
                first_guess = first_guess_1
                gauss_params = fit_gauss_2D(psf, first_guess = first_guess)
                # fitted center of the Gaussian
                psf_cx = gauss_params[1]
                psf_cy = gauss_params[2]
            except RuntimeError:
                print('WARNING: second Gaussian fit failed.')
                print('ERROR: Off-axis PSF could not be properly fitted.'+
                      '\n WARNING: input first guess position is used.')
                psf_cx = cyx_guess[1]
                psf_cy = cyx_guess[0]
        
        psf_coord[i] = (psf_cy, psf_cx)
    
    return psf_coord

def get_aperture_flux(img, radius_pix, cx=None, cy=None, t_int=1.):
    """ 
    Estimates the flux in an aperture of the input image by integrating the 
    pixel counts in a circular area centered on cx,cy and of radius radiux_pix.
    The circular area is defined by a whole number of pixels (no pixel fraction
    taken into account).
    Integration time of the image can also be provided in order to normalize 
    the final flux estimate.

    Parameters
    ----------
    psf : array_like
        Input 2D array.
    radius_pix : float
        Radius of the mask in pixels.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    t_int : float, optional
        Integration time of the psf image.

    Returns
    -------
    psf_flux : float
        Flux of the psf integrated in the circular area and normalized by its
        integration time if provided.
    """

    ny, nx = img.shape

    if cx == None :
        cx = (nx-1.) / 2.
    if cy == None :
        cy = (ny-1) / 2.

    circ_mask = circle_mask(img, radius_pix, cx=cx, cy=cy)
    flux = np.sum(img*circ_mask)
    
    return flux / t_int


def gauss_2D(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """ 
    Model function of a 2D Gaussian profile.
    """    
    x=xy[0]
    y=xy[1]
    
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def fit_gauss_2D(img, first_guess=None):
    ''' 
    Fits a 2D Gaussian pattern on the image.
    Returns the best fit parameters of the Gaussian shape.
    See model function gauss_2D(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''

    nx=img.shape[1]
    ny=img.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    x, y = np.meshgrid(x, y)

    init_xmax=np.unravel_index(img.argmax(), img.shape)[1]
    init_ymax=np.unravel_index(img.argmax(), img.shape)[0]
    
    if first_guess is None:
        first_guess = (img.max(), init_xmax, init_ymax, 5, 5, 0, 0)
        
    popt, pcov = curve_fit(gauss_2D, (x, y), img.ravel(), 
                               p0=first_guess)
    
    return popt

def display_tiptilt_target(img, image_sampling, vortex_cyx, 
                             tt_estimate_xy,
                             subim_width_lbdd=3., fig_num=1,
                             img_circ_rad=None, tt_circ_rad=None,
                             tt_lim = 1., 
                             plot_title=''):
    """
    Display the coronagraphic image and the corresponding tip-tilt estimate as
    a vector.
    """
    
    subim_width_pix  = subim_width_lbdd * image_sampling
    
    subim = subimage(img, vortex_cyx[1], vortex_cyx[0], subim_width_pix)
    
    subim_extent = np.array([-subim_width_lbdd, subim_width_lbdd,
                             -subim_width_lbdd, subim_width_lbdd])
    
    plt.figure(num=fig_num, figsize=(9,4))
    plt.clf()
    plt.suptitle(plot_title)
    
    fig, ax = plt.subplots(nrows=1, ncols=2, num=fig_num)
    
    ### LEFT PLOT: SUB-IMAGE
    ax[0].set_title('Coronagraphic image')
    ax[0].set_xlabel(r'$\lambda/D$')
    ax[0].set_ylabel(r'$\lambda/D$')
    ax[0].imshow(subim, cmap='plasma', 
               extent = subim_extent, 
               origin = 'lower', interpolation=None)
    if img_circ_rad is not None:
        for r in img_circ_rad:
            thet = np.linspace(0., 2*np.pi)
            ax[0].plot(r*np.cos(thet),r*np.sin(thet), 'w:')
    
    ### RIGHT PLOT: TIPTILT ESTIMATE
    ax[1].set_yticklabels([])
    ax1 = ax[1].twinx()
    ax1.vlines(0, -tt_lim, tt_lim, linestyle='--', color='.7')
    ax1.hlines(0, -tt_lim, tt_lim, linestyle='--', color='.7')
    if tt_circ_rad is not None:
        for r in tt_circ_rad:
            thet = np.linspace(0., 2*np.pi)
            ax1.plot(r*np.cos(thet),r*np.sin(thet), ':', color='.8')
    ax1.axis('equal')
    ax1.set_title('Tip-tilt estimate')
    ticks = np.linspace(-tt_lim, tt_lim, 5)
    ax1.set_xticks(ticks)
    ax1.set_yticks(ticks)
    ax1.plot([0., tt_estimate_xy[0]], [0., tt_estimate_xy[1]], 'r', linewidth=4)

    ax1.set_xlabel(r'$\lambda/D$')
    ax1.set_ylim(-tt_lim,tt_lim)
    ax1.set_xlim(-tt_lim,tt_lim)
    ax1.set_ylabel(r'$\lambda/D$')    

    return fig, ax

def display_tiptilt_sequence(tt_xy, delta_t = None, 
                             tt_circ_rad=None, fignum=2):
    """
    Display the tip-tilt estimates over time.
    """
    plt.figure(num=fignum, figsize=(9,4))
    plt.clf()
    plt.suptitle('Tip-tilt estimates over time')
    
    fig, ax = plt.subplots(nrows=1, ncols=2, num=fignum)
    
    ### LEFT PLOT: 2D tip-tilt estimates
    ax[0].plot(tt_xy[:,0],tt_xy[:,1], 'c', alpha=.3)
    ax[0].plot(tt_xy[:,0],tt_xy[:,1], 'co', alpha=.5, markersize=3)
    ax[0].axis('equal')
    ax[0].axis([-1.,1.,-1.,1.])
    ticks = np.linspace(-1., 1., 11)
    ax[0].set_xticks(ticks)
    ax[0].set_yticks(ticks)
    ax[0].vlines(0., -1., 1., color='.65', linestyle='--')
    ax[0].hlines(0., -1., 1., color='.65', linestyle='--')
    ax[0].set_xlabel(r'$\lambda/D$')
    ax[0].set_ylabel(r'$\lambda/D$')
    if tt_circ_rad is not None:
        for r in tt_circ_rad:
            thet = np.linspace(0., 2*np.pi)
            ax[0].plot(r*np.cos(thet),r*np.sin(thet), ':', color='.8')
    
    ### RIGHT PLOT: tip-tilt amplitude as a function of time ###
    n_tt      = tt_xy.shape[0]
    if delta_t is None:
        delta_t = 1.
        xlabel = '# frames'
    else:
        xlabel = 'Time'
    time_axis = np.linspace(0., n_tt*delta_t, n_tt)
    tt_mod    = np.sqrt(np.sum(tt_xy**2, axis=1))
    ymax      = np.ceil(np.max(tt_mod)*1.2*10.)/10.
    #yticks    = np.linspace(0, ymax, num=10)
    #ax[1].set_yticks(yticks)
    #ax[1].set_yticklabels([])
    ax[1].set_xlabel(xlabel)
    ax1 = ax[1].twinx()
    ax1.grid(color='.8', linestyle='--')
    #ax[1].set_ylim(0,ymax)
    ax1.plot(time_axis, tt_mod, 'co', alpha=.5, markersize=4)

    ax1.set_ylabel(r'Tip-tilt amplitude [$\lambda/D$]')
    ax1.set_xlim(0,)
    
    return fig, ax
    
def display_tiptilt_sequence_alma(tt_xy, delta_t = None, 
                             tt_circ_rad=None, tt_xy_simul = None, fignum=2):
    """
    Display the tip-tilt estimates over time.
    """
    plt.figure(num=fignum, figsize=(14.3,4)) #figsize=(9,4))
    plt.clf()
    if tt_xy_simul is None:
        plt.suptitle('Tip-tilt estimates over time')
    else:
        plt.suptitle('Tip-tilt measured-simulated differences over time')
    
    fig, ax = plt.subplots(nrows=1, ncols=3, num=fignum)
    
    ### LEFT PLOT: 2D tip-tilt estimates
    if tt_xy_simul is None:
        #ax[0].plot(tt_xy[:,0],tt_xy[:,1], 'c', color='blue', alpha=.3)
        ax[0].plot(tt_xy[:,0],tt_xy[:,1], 'co', color='blue', alpha=.5, markersize=3)
    else:
        #ax[0].plot(tt_xy_simul[:,0],tt_xy_simul[:,1], 'c', color='red', alpha=.3)
        #ax[0].plot(tt_xy_simul[:,0],tt_xy_simul[:,1], 'co', color='red', alpha=.5, markersize=3)
        ax[0].plot((tt_xy-tt_xy_simul)[:,0],(tt_xy-tt_xy_simul)[:,1], 'co', color='red', alpha=.5, markersize=3)
    
    if tt_xy_simul is not None:
        plotrange=[-1,1]
        plotrange=[-.2,.2]
        #plotrange=[-.1,.1] #zoom for CNRS plot
    else:
        plotrange=[np.min(tt_xy[:,0])-.1,np.max(tt_xy[:,1])+.1] #[-.2,.2]
        plotrange=[-.1,.1]
    #print('plot range: ',plotrange)
    ax[0].axis('equal')
    ax[0].axis([plotrange[0],plotrange[1],plotrange[0],plotrange[1]])
    ticks = np.linspace(plotrange[0], plotrange[1], 5) #11)
    ax[0].set_xticks(ticks)
    ax[0].set_yticks(ticks)
    ax[0].vlines(0., plotrange[0], plotrange[1], color='.65', linestyle='--')
    ax[0].hlines(0., plotrange[0], plotrange[1], color='.65', linestyle='--')
    ax[0].set_xlabel(r'$\lambda/D$')
    ax[0].set_ylabel(r'$\lambda/D$')
    if tt_circ_rad is not None:
        for r in tt_circ_rad:
            thet = np.linspace(0., 2*np.pi)
            ax[0].plot(r*np.cos(thet),r*np.sin(thet), ':', color='.8')
    
    ### MIDDLE PLOT: tip-tilt amplitude as a function of time ###
    n_tt      = tt_xy.shape[0]
    if delta_t is None:
        delta_t = 1.
        xlabel = '# frames'
    else:
        xlabel = 'Time'
    time_axis = np.linspace(0., n_tt*delta_t, n_tt)
    tt_mod    = np.sqrt(np.sum(tt_xy**2, axis=1))
    if tt_xy_simul is None:
        ymin=0.
        ymax      = int(np.ceil(np.ceil(np.max(tt_mod)*1.2*10.)/10.))
        #print('ymax',ymax)
        yticks    = np.linspace(0, ymax, int(np.ceil(ymax/0.2)))
    else:
        tt_mod_simul=np.sqrt(np.sum(tt_xy_simul**2, axis=1))
        #ymax      = np.ceil(np.max([tt_mod,tt_mod_simul])*1.1*10.)/10.
        ymin      = int(np.ceil(np.min([tt_mod-tt_mod_simul])*1.1*10.)/10.)
        ymax      = int(np.ceil(np.max([tt_mod-tt_mod_simul])*1.1*10.)/10.)
        #print('ymax-ymin',ymax-ymin)
        yticks    = np.linspace(ymin, ymax, int(np.ceil((ymax-ymin)/0.05)))
        #print('ymax-ymin',ymax-ymin)
    ax[1].yaxis.set_visible(False) #Remove y axes before remake them (otherwise misaligned ticks btw 2 axes)
    ax[1].set_yticks(yticks)
    ax[1].set_yticklabels([])
    ax[1].set_xlabel(xlabel)
    ax1 = ax[1].twinx()
    ax1.grid(color='.8', linestyle='--')
    if tt_xy_simul is None:
        ax1.plot(time_axis, tt_mod, 'co', color='blue', alpha=.5, markersize=4)
    else:
        #Plot simulated tip-tilts if relevant
        # ax1.plot(time_axis, tt_mod_simul, 'co', color='red', alpha=.5, markersize=4)
        ax1.plot(time_axis, tt_mod-tt_mod_simul, 'co', color='red', alpha=.5, markersize=4)

    ax1.set_ylabel(r'Tip-tilt amplitude [$\lambda/D$]')
    ax1.set_xlim(0,)
    #ylimmax=1.1
    #ax1.set_ylim(ymin,ymax)
    #ax1.text(7, ymax*.94, 'Measured',color='blue')
    #ax1.text(7, ymax*.89, 'Simulated',color='red')
    
    #ax1.tick_params(
    #axis='y',          # changes apply to the x-axis
    #which='both',      # both major and minor ticks are affected
    #left=False)      # ticks along the bottom edge are off
    
    ### RIGHT PLOT: tip-tilt amplitude as a function of time ###
    tt_dir=np.arctan2(tt_xy[:,1],tt_xy[:,0]) * 180./np.pi
    if tt_xy_simul is None:
        ymin=-180.
        ymax=180.
        yticks    = np.linspace(ymin, ymax, int(np.ceil(ymax/45.)))
    else:
        ymax=180. #60. #360.
        ymin=-180. #-60. #-360.
        yticks    = np.linspace(ymin, ymax, int(np.ceil((ymax-ymin)/45.))) #10. #45.
    ax[2].yaxis.set_visible(False) #Remove y axes before remake them (otherwise misaligned ticks btw 2 axes)
    ax[2].set_yticks(yticks)
    ax[2].set_yticklabels([])
    ax[2].set_xlabel(xlabel)
    ax2 = ax[2].twinx()
    ax2.grid(color='.8', linestyle='--')
    ax2.set_yticks(yticks)
    if tt_xy_simul is None:
        ax2.plot(time_axis, tt_dir, 'co', color='blue', alpha=.5, markersize=4)
    else:
        tt_dir_simul=np.arctan2(tt_xy_simul[:,1],tt_xy_simul[:,0]) * 180./np.pi
        #tt_dir_diff=(tt_dir-tt_dir_simul) % (np.sign(tt_dir-tt_dir_simul)*360)
        tt_dir_diff=((tt_dir-tt_dir_simul)+180.) % (np.sign(tt_dir-tt_dir_simul+180.)*360)-np.sign((tt_dir-tt_dir_simul)+180.)*180.
            
        #print('tt_dir simul: ',tt_dir_simul)
        #print('tt_dir meas: ',tt_dir)
        #print('tt_dir meas-simul: ',tt_dir_diff)
        #Plot simulated tip-tilts if relevant
        # ax2.plot(time_axis, np.arctan2(tt_xy_simul[:,1],tt_xy_simul[:,0]) * 180./np.pi, 'co', color='red', alpha=.5, markersize=4)
        ax2.plot(time_axis, tt_dir_diff, 'co', color='red', alpha=.5, markersize=4)


    ax2.set_ylabel(r'Tip-tilt direction [deg]')
    ax2.set_xlim(0.,)
    ax2.set_ylim(ymin,ymax)
    #ax2.text(7, ymax*.94, 'Measured',color='blue')
    #ax2.text(7, ymax*.89, 'Simulated',color='red')
    
    #plt.subplots_adjust(wspace=0.3)
    plt.tight_layout()
    fig.savefig("Figure_42.pdf")

    return fig, ax

def display_tiptilt_sequence_alma2(tt_xy, delta_t = None,
                             tt_circ_rad=None, tt_xy_simul = None, fignum=2):
    """
    Display the tip-tilt estimates over time.
    """
    #plt.figure(num=fignum, figsize=(14.3/3.,4)) #figsize=(9,4))
    plt.clf()
    #fig, ax = plt.subplots(num=fignum, figsize = (7,7),dpi=200)
    #fig, ax = plt.subplots(num=fignum, figsize = (6,6),dpi=200)
    fig, ax = plt.subplots(num=fignum, figsize = (7,7),dpi=200)
    #plt.rcParams["font.size"] = 18
    #plt.rcParams.update({'font.size': 13})
    #plt.rc('xtick', labelsize=13) 
    #plt.rc('ytick', labelsize=13)
    #if tt_xy_simul is None:
    #    plt.suptitle('Tip-tilt estimates over time')
    #else:
    #    plt.suptitle('Tip-tilt measured-simulated differences over time')

    #fig, ax = plt.subplots(nrows=1, ncols=1, num=fignum)

    ### PLOT: 2D tip-tilt estimates
    if tt_xy_simul is None:
        #ax[0].plot(tt_xy[:,0],tt_xy[:,1], 'c', color='blue', alpha=.3)
        plt.plot(tt_xy[:,0],tt_xy[:,1], 'co', color='blue', alpha=.5, markersize=3)
    else:
        #ax[0].plot(tt_xy_simul[:,0],tt_xy_simul[:,1], 'c', color='red', alpha=.3)
        #ax[0].plot(tt_xy_simul[:,0],tt_xy_simul[:,1], 'co', color='red', alpha=.5, markersize=3)
        plt.plot((tt_xy-tt_xy_simul)[:,0],(tt_xy-tt_xy_simul)[:,1], 'co', color='red', alpha=.5, markersize=3)

    if tt_xy_simul is not None:
        plotrange=[-1,1]
        plotrange=[-.2,.2]
        #plotrange=[-.1,.1] #zoom for CNRS plot
    else:
        plotrange=[np.min(tt_xy[:,0])-.1,np.max(tt_xy[:,1])+.1] #[-.2,.2]
        plotrange=[-.1,.1]
    #print('plot range: ',plotrange)
    #plt.axis('equal')
    plt.axis([plotrange[0],plotrange[1],plotrange[0],plotrange[1]])
    ticks = np.linspace(plotrange[0], plotrange[1], 5) #11)
    print('ticks',ticks)
    #ax.set_xticklabels(ticks, fontsize=12)
    #ax.set_yticklabels(ticks, fontsize=12)
    ftsz=16 #11 #13
    plt.xticks(fontsize=ftsz)
    plt.yticks(fontsize=ftsz)
    #plt.setp(ax.get_xticklabels(), fontsize=12)
    #plt.setp(ax.get_yticklabels(), fontsize=12)
    ax.set_xticks(ticks) #,labelsize=13)
    ax.set_yticks(ticks) #,labelsize=13)
    #ax.set_xticklabels(ticks, fontsize=12)
    #ax.set_yticklabels(ticks, fontsize=12)
    plt.xlabel(r'$\lambda/D$',fontsize=ftsz)
    plt.ylabel(r'$\lambda/D$',fontsize=ftsz)
    plt.vlines(0., plotrange[0], plotrange[1], color='.65', linestyle='--')
    plt.hlines(0., plotrange[0], plotrange[1], color='.65', linestyle='--')
    ax.set_xlabel(r'$\lambda/D$')
    ax.set_ylabel(r'$\lambda/D$')
    if tt_circ_rad is not None:
        for r in tt_circ_rad:
            thet = np.linspace(0., 2*np.pi)
            plt.plot(r*np.cos(thet),r*np.sin(thet), ':', color='.8')

    #plt.subplots_adjust(wspace=0.9)
    plt.tight_layout()
    fig.savefig("Figure_42.pdf")

    return fig, ax

def dist(NAXIS):
    """Returns a rectangular array in which the value of each element is proportional to its frequency.
    Author: Abeelen
    >>> dist(3)
    array([[ 0.        ,  1.        ,  1.        ],
           [ 1.        ,  1.41421356,  1.41421356],
           [ 1.        ,  1.41421356,  1.41421356]])
    >>> dist(4)
    array([[ 0.        ,  1.        ,  2.        ,  1.        ],
           [ 1.        ,  1.41421356,  2.23606798,  1.41421356],
           [ 2.        ,  2.23606798,  2.82842712,  2.23606798],
           [ 1.        ,  1.41421356,  2.23606798,  1.41421356]])
    """
    axis = np.linspace(-NAXIS/2+1, NAXIS/2, NAXIS)
    result = np.sqrt(axis**2 + axis[:,np.newaxis]**2)
    return np.roll(result, int(NAXIS/2+1), axis=(0,1))
