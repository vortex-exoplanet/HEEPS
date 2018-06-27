#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:01:50 2018

@author: brunella
"""

import proper
import numpy as np
import math
from astropy.io import fits
import cv2
from circular_apodization import circular_apodization

def simple_telescope(wavelength, gridsize,PASSVALUE = {'prefix':'prova', 'charge':0, 'diam':37., 'beam_ratio':0.25,  'RAVC':True, 'LS':True}):
    
     ## call all the vues passed via passvalue
    prefix = PASSVALUE['prefix']
    charge = PASSVALUE['charge']
    diam = PASSVALUE['diam']
    beam_ratio = PASSVALUE['beam_ratio']
    RAVC = PASSVALUE['RAVC']
    LS = PASSVALUE['LS']
    

    n = int(gridsize)
    npupil = math.ceil(n*beam_ratio) # compute the pupil size --> has to be ODD (proper puts the center in the up right pixel next to the grid center)
    if npupil % 2 == 0:
        npupil = npupil +1
    print("npupil: ", npupil)
    r_obstr = 0.3
    f_lens = 17.8 * diam
     
        
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio) # define the simualtion pupil
    lamda=proper.prop_get_wavelength(wfo) #save the wavelength value [m] into lamda
    
    print("lambda: ", lamda)
    
    proper.prop_circular_aperture(wfo, diam/2)
    proper.prop_circular_obscuration(wfo, r_obstr, NORM=True)
    proper.prop_rectangular_obscuration(wfo, 0.60, 2*diam,ROTATION=0.) # define the spiders
    proper.prop_rectangular_obscuration(wfo, 0.60, 2*diam,ROTATION=60.) # define the spiders
    proper.prop_rectangular_obscuration(wfo, 0.60, 2*diam,ROTATION=120.) # define the spiders

    fits.writeto(prefix+'_pupil_pre_define.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)

    proper.prop_define_entrance(wfo) #define the entrance wavefront
    
    print("pupil ")

    
    if (RAVC == True): # when tha apodizer is present
        t1_opt = 1. - 1./4*(r_obstr**2 + r_obstr*(math.sqrt(r_obstr**2 + 8.))) # define the apodizer transmission [Mawet2013]
        R1_opt = (r_obstr/math.sqrt(1. - t1_opt)) # define the apodizer radius [Mawet2013]
        print ("r1_opt: ", R1_opt)
        print ("t1_opt: ", t1_opt)
        apodizer=circular_apodization(wfo, R1_opt, 1., t1_opt, NORM=True) # define the apodizer
        apodizer = proper.prop_shift_center(apodizer)
        proper.prop_multiply(wfo, apodizer)
    
    print("RAVC ")


    fits.writeto(prefix+'_pupil_amplitude_charge'+str(charge)+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
    fits.writeto(prefix+'_pupil_phase_charge'+str(charge)+'.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
 
    proper.prop_propagate(wfo, f_lens, 'inizio') # propagate wavefront

    proper.prop_lens(wfo, f_lens, 'focusing lens vortex') # propagate through a lens
    proper.prop_propagate(wfo, f_lens, 'VC') # propagate wavefront

    print("pre-vortex ")

    ofst = 0 # no offset
    ramp_sign = 1 #sign of charge is positive
    sampling = n
    ramp_oversamp = 11. # vortex is oversampled for a better discretization
    nramp = int(n*ramp_oversamp) #oversamp
            # create the vortex by creating a matrix (theta) representing the ramp (created by atan 2 gradually varying matrix, x and y)
    y1 = np.ones((nramp,), dtype=np.int)
    y2 = np.arange(0, nramp, 1.) - (nramp/2) - int(ramp_oversamp)/2
    y = np.outer(y2, y1)
    x = np.transpose(y)
    theta = np.arctan2(y,x)
    x = 0
    y = 0
    vvc_tmp = np.exp(1j*(ofst + ramp_sign*charge*theta))
    theta = 0
    vvc_real_resampled = cv2.resize(vvc_tmp.real, (0,0), fx=1/ramp_oversamp, fy=1/ramp_oversamp, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
    vvc_imag_resampled = cv2.resize(vvc_tmp.imag, (0,0), fx=1/ramp_oversamp, fy=1/ramp_oversamp, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
    vvc = np.array(vvc_real_resampled, dtype=complex)
    vvc.imag = vvc_imag_resampled
    vvcphase = np.arctan2(vvc.imag, vvc.real) # create the vortex phase
    vvc_complex = np.array(np.zeros((n,n)), dtype=complex)
    vvc_complex.imag = vvcphase
    vvc = np.exp(vvc_complex)
    vvc = proper.prop_shift_center(vvc)
    wfo._wfarr = wfo._wfarr*vvc # the wavefront takes into account the real pupil with the perfect-result vortex field

    print("post-vortex ")

    fits.writeto(prefix+'_afterVortex.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
    fits.writeto(prefix+'_afterVortex.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)

    proper.prop_propagate(wfo, f_lens, 'Lyot Collimetor') # propagate wavefront

    proper.prop_lens(wfo, f_lens, 'Lyot Collimetor') # propagate wavefront through  a lens
    proper.prop_propagate(wfo, f_lens, 'Lyot Stop') # propagate wavefront


    fits.writeto(prefix+'_before_LS_charge'+str(charge)+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
  
    if (RAVC==True):
        t1_opt = 1. - 1./4*(r_obstr**2 + r_obstr*(math.sqrt(r_obstr**2 + 8.))) # define the apodizer transmission [Mawet2013]
        R1_opt = (r_obstr/math.sqrt(1. - t1_opt)) # define the apodizer radius [Mawet2013]
        r_LS = R1_opt + 0.03
    else:
        r_LS = r_obstr + 0.03

    print("pre-LS ")

    proper.prop_circular_aperture(wfo, 0.98, NORM=True)
    proper.prop_circular_obscuration(wfo, r_LS, NORM=True)
    proper.prop_rectangular_obscuration(wfo, 1.10, 2*diam,ROTATION=0.) # define the spiders
    proper.prop_rectangular_obscuration(wfo, 1.10, 2*diam,ROTATION=60.) # define the spiders
    proper.prop_rectangular_obscuration(wfo, 1.10, 2*diam,ROTATION=120.) # define the spiders

    print("post-LS ")


    fits.writeto(prefix+'_after_LS'+str(int(LS))+'_charge'+str(charge)+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)


    proper.prop_propagate(wfo, f_lens) # propagate wavefront

    proper.prop_lens(wfo, f_lens) # propagate wavefront through a lens
    proper.prop_propagate(wfo, f_lens) # propagate wavefront
    
    (wfo, sampling) = proper.prop_end(wfo, NOABS = True) # conclude the simulation --> noabs= the wavefront array will be complex
    
    
    return (wfo, sampling) # return the wavefront




