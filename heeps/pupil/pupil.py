import numpy as np
import cv2
import proper
from astropy.io import fits


def pupil(conf):

    diam = conf['DIAM']
    gridsize = conf['GRIDSIZE']
    spiders_width = conf['SPIDERS_WIDTH']
    spiders_angle = conf['SPIDERS_ANGLE']
    pixelsize = conf['PIXEL_SCALE'] 
    
    r_obstr = conf['R_OBSTR']
    wavelength = conf['WAVELENGTH']
    missing_segments_number = conf['MIS_SEGMENTS_NU']
    
    Debug = conf['DEBUG']
    Debug_print = conf['DEBUG_PRINT'] 

    prefix = conf['PREFIX']
    input_dir = conf['INPUT_DIR']
    out_dir = conf['OUT_DIR'] 
    
    beam_ratio = pixelsize*4.85e-9/(wavelength/diam)
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)
    n = int(gridsize)
    npupil = np.ceil(gridsize*beam_ratio) # compute the pupil size --> has to be ODD (proper puts the center in the up right pixel next to the grid center)
    if npupil % 2 == 0:
        npupil = npupil +1
    conf['NPUPIL'] = npupil
    if (Debug_print == True):
        print ("npupil: ", npupil)
        print("lambda: ", wavelength)
    pupil_file = fits.getdata(input_dir + conf['PUPIL_FILE'])
    if (missing_segments_number == 0):
        if (isinstance(pupil_file, (list, tuple, np.ndarray)) == True):
            pupil = pupil_file
            pupil_pixels = (pupil.shape)[0]## fits file size
            scaling_factor = float(npupil)/float(pupil_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
            if (Debug_print==True):
                print ("scaling_factor: ", scaling_factor)
            pupil_scale = cv2.resize(pupil.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
            if (Debug_print==True):
                print ("pupil_resample", pupil_scale.shape)
            pupil_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
            if (Debug_print==True):
                print("n: ", n)
                print("npupil: ", npupil)
            pupil_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] = pupil_scale # insert the scaled pupil into the 0s grid

        proper.prop_circular_aperture(wfo, diam/2) # create a wavefront with a circular pupil
    
        if (isinstance(pupil_file, (list, tuple, np.ndarray)) == True):
            proper.prop_multiply(wfo, pupil_large) # multiply the saved pupil
        else:
            proper.prop_circular_obscuration(wfo, r_obstr, NORM=True) # create a wavefront with a circular central obscuration
        if (spiders_width != 0):
            for iter in range(0,len(spiders_angle)):
                proper.prop_rectangular_obscuration(wfo, spiders_width, 2*diam, ROTATION = spiders_angle[iter]) # define the spiders    
    else:
        if (missing_segments_number == 1):
            pupil = fits.getdata(input_dir+'/ELT_2048_37m_11m_5mas_nospiders_1missing_cut.fits')
        if (missing_segments_number == 2):
            pupil = fits.getdata(input_dir+'/ELT_2048_37m_11m_5mas_nospiders_2missing_cut.fits')
        if (missing_segments_number == 4):
            pupil = fits.getdata(input_dir+'/ELT_2048_37m_11m_5mas_nospiders_4missing_cut.fits')
        if (missing_segments_number == 7):
            pupil = fits.getdata(input_dir+'/ELT_2048_37m_11m_5mas_nospiders_7missing_1_cut.fits')

        pupil_pixels = (pupil.shape)[0]## fits file size
        scaling_factor = float(npupil)/float(pupil_pixels) ## scaling factor between the fits file size and the pupil size of the simulation
        if (Debug_print==True):
            print ("scaling_factor: ", scaling_factor)
        pupil_scale = cv2.resize(pupil.astype(np.float32), (0,0), fx=scaling_factor, fy=scaling_factor, interpolation=cv2.INTER_LINEAR) # scale the pupil to the pupil size of the simualtions
        if (Debug_print==True):
            print ("pupil_resample", pupil_scale.shape)
        pupil_large = np.zeros((n,n)) # define an array of n-0s, where to insert the pupuil
        if (Debug_print==True):
            print("n: ", n)
            print("npupil: ", npupil)
        pupil_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] =pupil_scale # insert the scaled pupil into the 0s grid

            
        proper.prop_multiply(wfo, pupil_large) # multiply the saved pupil
        if (spiders_width!=0):
            for iter in range(0,len(spiders_angle)):
                proper.prop_rectangular_obscuration(wfo, spiders_width, 2*diam,ROTATION=spiders_angle[iter]) # define the spiders

    if (Debug==True):
        fits.writeto(out_dir + prefix +'_intial_pupil.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        
    proper.prop_define_entrance(wfo) #define the entrance wavefront
    wfo.wfarr *= 1./np.amax(wfo._wfarr) # max(amplitude)=1
    return wfo





