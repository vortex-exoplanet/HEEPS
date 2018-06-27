import numpy as np
import cv2
import proper
import math
from astropy.io import fits

def atmosphere(wf, npupil, atm_screen, Debug_print, Debug):
    
    n = int(proper.prop_get_gridsize(wf))
    lamda=proper.prop_get_wavelength(wf) #save the wavelength value [m] into lamda

    screen = atm_screen # read the phase screen
    if (Debug_print == True):
        print ("screen", screen.shape)
    screen_pixels = (screen.shape)[0] # size of the phase screen
    scaling_factor_screen = float(npupil)/float(screen_pixels) # scaling factor the phase screen to the simulation pupil size
    if (Debug_print == True):
        print ("scaling_factor_screen: ", scaling_factor_screen)
    screen_scale = cv2.resize(screen.astype(np.float32), (0,0), fx=scaling_factor_screen, fy=scaling_factor_screen, interpolation=cv2.INTER_LINEAR) # scale the the phase screen to the simulation pupil size
    if (Debug_print == True):
        print ("screen_resample", screen_scale.shape)
    screen_large = np.zeros((n,n)) # define an array of n-0s, where to insert the screen
    if (Debug_print == True):
        print("n: ", n)
        print("npupil: ", npupil)
    screen_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] =screen_scale # insert the scaled screen into the 0s grid


    lambda2=lamda/(1e-6) # need to use lambda in microns
    screen_atm = np.exp(1j*screen_large/lambda2*2*math.pi)
    proper.prop_multiply(wf, screen_atm) # multiply the atm screen to teh wavefront


    return wf
