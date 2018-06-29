import numpy as np
import proper
import math
import cv2
from circular_apodization import circular_apodization




def apodization(wf, r_obstr, npupil, RAVC=False, phase_apodizer_file=0, amplitude_apodizer_file=0, apodizer_misalignment=0, Debug_print=False):
    
    n = int(proper.prop_get_gridsize(wf))
    apodizer = 1
    if (RAVC == True):
        t1_opt = 1. - 1./4*(r_obstr**2 + r_obstr*(math.sqrt(r_obstr**2 + 8.))) # define the apodizer transmission [Mawet2013]
        R1_opt = (r_obstr/math.sqrt(1. - t1_opt)) # define the apodizer radius [Mawet2013]
        if (Debug_print == True):
            print ("r1_opt: ", R1_opt)
            print ("t1_opt: ", t1_opt)
        apodizer=circular_apodization(wf, R1_opt, 1., t1_opt, xc = apodizer_misalignment[0], yc = apodizer_misalignment[1],NORM=True) # define the apodizer
        apodizer = proper.prop_shift_center(apodizer)

    if (isinstance(phase_apodizer_file, (list, tuple, np.ndarray)) == True):
        xc_pixels = int(apodizer_misalignment[3]*npupil)
        yc_pixels = int(apodizer_misalignment[4]*npupil)
        apodizer_pixels = (phase_apodizer_file.shape)[0]## fits file size
        scaling_factor = float(npupil)/float(apodizer_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
        if (Debug_print==True):
            print ("scaling_factor: ", scaling_factor)
        apodizer_scale = cv2.resize(phase_apodizer_file.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
        if (Debug_print==True):
            print ("apodizer_resample", apodizer_scale.shape)
        apodizer_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
        if (Debug_print==True):
            print("n: ", n)
            print("npupil: ", npupil)
        apodizer_large[int(n/2)+1-int(npupil/2)-1 + xc_pixels:int(n/2)+1+int(npupil/2)+ xc_pixels,int(n/2)+1-int(npupil/2)-1+ yc_pixels:int(n/2)+1+int(npupil/2)+ yc_pixels] =apodizer_scale # insert the scaled pupil into the 0s grid
        apodizer = np.exp(1j*apodizer_large)


    if (isinstance(amplitude_apodizer_file, (list, tuple, np.ndarray)) == True):
        xc_pixels = int(apodizer_misalignment[0]*npupil)
        yc_pixels = int(apodizer_misalignment[1]*npupil)
        apodizer_pixels = (amplitude_apodizer_file.shape)[0]## fits file size
        scaling_factor = float(npupil)/float(apodizer_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
        if (Debug_print==True):
            print ("scaling_factor: ", scaling_factor)
        apodizer_scale = cv2.resize(amplitude_apodizer_file.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
        if (Debug_print==True):
            print ("apodizer_resample", apodizer_scale.shape)
        apodizer_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
        if (Debug_print==True):
            print("n: ", n)
            print("npupil: ", npupil)
        apodizer_large[int(n/2)+1-int(npupil/2)-1 + xc_pixels:int(n/2)+1+int(npupil/2)+ xc_pixels,int(n/2)+1-int(npupil/2)-1+ yc_pixels:int(n/2)+1+int(npupil/2)+ yc_pixels] =apodizer_scale # insert the scaled pupil into the 0s grid
        apodizer = apodizer_large
    proper.prop_multiply(wf, apodizer)
    
    return wf


