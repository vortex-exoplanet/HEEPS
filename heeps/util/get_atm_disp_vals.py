from astropy.io import fits 
import numpy as np

def get_disp_vals(adc_nu, steps, **conf):
    ''' Function to load the pre-calculated dispersion values
        and interpolate it to give dispersion values at any zenith angle
    Input: ADC number, wavelength steps to calculate the dispersion values
       
    Output: wavelengths at which the calculation is done, fit for residual dispersion in 
            x-axis, fit for residual dispersion in y-axis, and atmospheric disp in y-axis.
    Note: the atmospheric dispersion in x-axis is 0.        
    
    '''
    # These are latest calcultion of disp done by Roy.
    # Optimized for center of the HCI-Llong band

    filename = 'PSF_offsets_HCI-Llong_method_scale_zOpt{}.0.fits'.format(adc_nu)
    hdul = fits.open(conf['dir_input']+'dispersion_files/'+filename)
    wavelength = hdul[0].data  # Sampled into 256 points
    zenith_angles = np.rad2deg(hdul[1].data) # Starts from 0 to 60 degrees (61 points)
    xoffsets = hdul[2].data # residual x-disp, in the presence of an ADC
    yoffsets = hdul[3].data # residual y-disp, in the presence of an ADC
    xdisp = hdul[4].data # atmospheric disp in x-direction (0)
    ydisp = hdul[5].data # atmospheric disp in y-direction (elevation axis)

    # Selecting HCI-L long band
    ind1, = np.where(wavelength==3.7)[0]
    ind2, = np.where(wavelength==3.9541178)[0]
#    L_hci_short = wavelength[132:170]  # 3.50 - 3.70µm
    L_hci_long = wavelength[ind1:ind2+1]  # 3.70 - 3.95µm
    ss = np.linspace(0, L_hci_long.shape[0]-1, steps, dtype=int) 
    sel_L_long = L_hci_long[ss]

    fit_order = 5 # this fits dispersion curves well (other fits were tried as well) 
    fit_res_x = np.zeros((steps, fit_order+1))
    fit_res_y = np.zeros((steps, fit_order+1))
    # this calculates the fitting parameters of resdiual disp for selected wavelengths
    for i, index in enumerate(ss+ind1):
        fit_res_x[i] = np.polyfit(zenith_angles, xoffsets[index,:], fit_order)
        fit_res_y[i] = np.polyfit(zenith_angles, yoffsets[index,:], fit_order)
    # Since atmospheric disp in x-axis is zero, only in y-axis is calculated
    fit_disp_y = np.zeros((steps, fit_order+1))
    for i, index in enumerate(ss+ind1):
        fit_disp_y[i] = np.polyfit(zenith_angles, ydisp[index,:], fit_order)
    return sel_L_long, fit_res_x, fit_res_y, fit_disp_y

def get_transit_para_1hr_case(**conf):
    '''Get the transit parameters '''
    data = np.loadtxt(conf['dir_input']+'dispersion_files/'+'ha_parang_zangle.txt')
    ha = data[:,0]
    para = data[:,1]
    star_zenith = data[:,2]
    return star_zenith

def get_transit_para_alpha_cen(**conf):
    '''Get the transit parameters '''
    data = np.loadtxt(conf['dir_input']+'dispersion_files/'+'ha_parang_zangle_alfCen.txt') # 6hr Alpha-case
    ha = data[:,0]
    para = data[:,1]
    star_zenith = data[:,2]
    return star_zenith
