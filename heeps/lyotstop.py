import numpy as np
import proper
import math
import cv2
from astropy.io import fits

def lyotstop(wf, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug):

    if (RAVC==True): # define the inner radius of the Lyot Stop
        t1_opt = 1. - 1./4*(r_obstr**2 + r_obstr*(math.sqrt(r_obstr**2 + 8.))) # define the apodizer transmission [Mawet2013]
        R1_opt = (r_obstr/math.sqrt(1. - t1_opt)) # define teh apodizer radius [Mawet2013]
        r_LS = R1_opt + LS_parameters[1] # when a Ring apodizer is present, the inner LS has to have at least the value of the apodizer radius
    else:
        r_LS = r_obstr + LS_parameters[1] # when no apodizer, the LS has to have at least the radius of the pupil central obstruction
    if LS==True: # apply the LS
        if (Debug_print==True):
            print("LS parameters: ", LS_parameters)
        proper.prop_circular_aperture(wf, LS_parameters[0], LS_misalignment[0], LS_misalignment[1], NORM=True)
        proper.prop_circular_obscuration(wf, r_LS, LS_misalignment[0], LS_misalignment[1], NORM=True)
        if (LS_parameters[2]!=0):
            for iter in range(0,len(spiders_angle)):
                if (Debug_print==True):
                    print("LS_misalignment: ", LS_misalignment)
                proper.prop_rectangular_obscuration(wf, LS_parameters[2], 2*diam,LS_misalignment[0], LS_misalignment[1], ROTATION=spiders_angle[iter]) # define the spiders
        if (Debug==True):
            n = proper.prop_get_gridsize(wf)
            out_dir = str('./output_files/')
            fits.writeto(out_dir +'_Lyot_stop.fits', proper.prop_get_amplitude(wf)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        
    if (isinstance(LS_phase_apodizer_file, (list, tuple, np.ndarray)) == True):
        xc_pixels = int(LS_misalignment[3]*npupil)
        yc_pixels = int(LS_misalignment[4]*npupil)
        apodizer_pixels = (LS_phase_apodizer_file.shape)[0]## fits file size
        scaling_factor = float(npupil)/float(pupil_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
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
        phase_multiply = np.array(np.zeros((n,n)), dtype=complex) # create a complex array
        phase_multiply.imag = apodizer_large # define the imaginary part of the complex array as the atm screen
        apodizer = np.exp(phase_multiply)
        proper.prop_multiply(wf, apodizer)
        if (Debug == True):
            fits.writeto('LS_apodizer.fits', proper.prop_get_phase(wf), overwrite=True)


    
    if (isinstance(LS_amplitude_apodizer_file, (list, tuple, np.ndarray)) == True):
        print('4th')
        xc_pixels = int(LS_misalignment[0]*npupil)
        yc_pixels = int(LS_misalignment[1]*npupil)
        apodizer_pixels = (LS_amplitude_apodizer_file.shape)[0]## fits file size
        scaling_factor = float(npupil)/float(pupil_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
        if (Debug_print==True):
            print ("scaling_factor: ", scaling_factor)
            apodizer_scale = cv2.resize(amplitude_apodizer_file.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
        if (Debug_print==True):
            print ("apodizer_resample", apodizer_scale.shape)
        apodizer_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
        if (Debug_print==True):
            print("grid_size: ", n)
            print("npupil: ", npupil)
        apodizer_large[int(n/2)+1-int(npupil/2)-1 + xc_pixels:int(n/2)+1+int(npupil/2)+ xc_pixels,int(n/2)+1-int(npupil/2)-1+ yc_pixels:int(n/2)+1+int(npupil/2)+ yc_pixels] =apodizer_scale # insert the scaled pupil into the 0s grid
        apodizer = apodizer_large
        proper.prop_multiply(wf, apodizer)
        if (Debug == True):
            fits.writeto('LS_apodizer.fits', proper.prop_get_amplitude(wf), overwrite=True)

    return wf
