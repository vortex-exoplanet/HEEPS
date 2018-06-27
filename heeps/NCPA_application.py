
import numpy as np
import cv2
import proper
import math
from astropy.io import fits

def NCPA_application(wf, npupil, NCPA, path, Debug_print, Debug):
    
    n = int(proper.prop_get_gridsize(wf))
    lamda=proper.prop_get_wavelength(wf) #save the wavelength value [m] into lamda
    
    if (Debug_print == True):
        print ("NCPA", NCPA.shape)
    NCPA_pixels = (NCPA.shape)[0] # size of the phase screen
    scaling_factor_NCPA = float(npupil)/float(NCPA_pixels) # scaling factor the phase screen to the simulation pupil size
    if (Debug_print == True):
        print ("scaling_factor_NCPA: ", scaling_factor_NCPA)
    NCPA_scale = cv2.resize(NCPA.astype(np.float32), (0,0), fx=scaling_factor_NCPA, fy=scaling_factor_NCPA, interpolation=cv2.INTER_LINEAR) # scale the the phase screen to the simulation pupil size
    if (Debug_print == True):
        print ("NCPA_resample", NCPA_scale.shape)
    NCPA_large = np.zeros((n,n)) # define an array of n-0s, where to insert the screen
    if (Debug_print == True):
        print("n: ", n)
        print("npupil: ", npupil)
    NCPA_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] =NCPA_scale # insert the scaled screen into the 0s grid
    if (Debug == 1):
        fits.writeto(path + 'NCPA_large.fits', NCPA_large, overwrite=True) # fits file for the screen
    
    
    lambda2=lamda/(1e-6) # need to use lambda in microns
    screen_atm = np.exp(1j*NCPA_large/lambda2*2*math.pi)
    proper.prop_multiply(wf, screen_atm) # multiply the atm screen to teh wavefront

#rms_error = NCPA[0]#30nm-- RMS wavefront error in meters
#c_freq = NCPA[1]# 5  ;-- correlation frequency (cycles/meter)
#high_power = NCPA[2]#3.0  ;-- high frequency falloff (r^-high_power)
#RANDOM_MAP = readfits('RANDOM_MAP.fits')
#proper.prop_psd_errormap_mod( wf, rms_error, c_freq, high_power, RMS=True)
    
    
    return
