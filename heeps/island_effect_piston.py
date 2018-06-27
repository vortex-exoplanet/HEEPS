import numpy as np
import cv2
import proper
import vip_hci
import math
from astropy.io import fits
import os

def island_effect_piston(wf, npupil, Island_Piston, path, Debug_print, Debug):

    n = int(proper.prop_get_gridsize(wf))
    lamda=proper.prop_get_wavelength(wf) #save the wavelength value [m] into lamda

    PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))
    petal1 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal1_243px.fits')
    petal2 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal2_243px.fits')
    petal3 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal3_243px.fits')
    petal4 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal4_243px.fits')
    petal5 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal5_243px.fits')
    petal6 = fits.getdata(PACKAGE_PATH+'/1024_pixelsize5mas_petal6_243px.fits')

    piston_petal1 = Island_Piston[0]*petal1
    piston_petal2 = Island_Piston[1]*petal2
    piston_petal3 = Island_Piston[2]*petal3
    piston_petal4 = Island_Piston[3]*petal4
    piston_petal5 = Island_Piston[4]*petal5
    piston_petal6 = Island_Piston[5]*petal6

    piston = piston_petal1 + piston_petal2 + piston_petal3 + piston_petal4 + piston_petal5 + piston_petal6

    piston_pixels = (piston.shape)[0]## fits file size
    scaling_factor = float(npupil)/float(piston_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
    if (Debug_print==True):
        print ("scaling_factor: ", scaling_factor)
    piston_scale = cv2.resize(piston.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
    if (Debug_print==True):
        print ("piston_resample", piston_scale.shape)
    piston_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
    if (Debug_print==True):
        print("n: ", n)
        print("npupil: ", npupil)
    piston_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] =piston_scale # insert the scaled pupil into the 0s grid
    if (Debug == 1):
        fits.writeto(path+'piston_phase.fits', piston_large, overwrite=True) # fits file for the screen


    lambda2=lamda/(1e-6) # need to use lambda in microns
    piston_phase = np.exp(1j*piston_large/lambda2*2*math.pi)
    proper.prop_multiply(wf, piston_phase) # multiply the atm screen to teh wavefront

    return
