import proper
import numpy as np
import os

from astropy.io import fits


def telescope(wavelength, gridsize,PASSVALUE = {'prefix':'prova', 'path':os.path.abspath(os.path.join(__file__, os.pardir)), 'charge':0, 'CAL':0, 'diam':37., 
                                                'spiders_width':0.60, 'spiders_angle':[0., 60., 120.], 'beam_ratio': 0.25, 'f_lens':658.6, 'npupil':243, 
                                                'r_obstr':0.3, 'pupil_file':0, 'phase_apodizer_file':0, 'amplitude_apodizer_file':0, 'TILT':[0., 0.],
                                                'LS':False,'RAVC':False, 'LS_phase_apodizer_file':0, 'LS_amplitude_apodizer_file':0,'LS_parameters':[0.0, 0.0, 0.0], 'atm_screen':0, 'missing_segments_number':0, 'apodizer_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                                'LS_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Island_Piston':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],'NCPA':0, 
                                                'Debug_print':False,'Debug':False}):
    
     ## call all the vues passed via passvalue
    prefix = PASSVALUE['prefix']
    path = PASSVALUE['path']
    charge = PASSVALUE['charge']
    CAL = PASSVALUE['CAL']
    diam = PASSVALUE['diam']
    spiders_width = PASSVALUE['spiders_width']
    spiders_angle = PASSVALUE['spiders_angle']
    beam_ratio = PASSVALUE['beam_ratio']
    f_lens = PASSVALUE['f_lens']
    npupil = PASSVALUE['npupil']
    r_obstr = PASSVALUE['r_obstr']
    pupil_file = PASSVALUE['pupil_file']
    phase_apodizer_file = PASSVALUE['phase_apodizer_file']
    amplitude_apodizer_file = PASSVALUE['amplitude_apodizer_file']
    TILT = PASSVALUE['TILT']
    LS = PASSVALUE['LS']
    RAVC = PASSVALUE['RAVC']
    LS_phase_apodizer_file = PASSVALUE['LS_phase_apodizer_file']
    LS_amplitude_apodizer_file = PASSVALUE['LS_amplitude_apodizer_file']
    LS_parameters = PASSVALUE['LS_parameters']
    atm_screen = PASSVALUE['atm_screen']
    missing_segments_number = PASSVALUE['missing_segments_number']
    apodizer_misalignment = PASSVALUE['apodizer_misalignment']
    LS_misalignment = PASSVALUE['LS_misalignment']
    Island_Piston = PASSVALUE['Island_Piston']
    NCPA = PASSVALUE['NCPA']
    Debug_print = PASSVALUE['Debug_print']
    Debug = PASSVALUE['Debug']
    
    TILT=np.array(TILT)
    apodizer_misalignment=np.array(apodizer_misalignment)
    LS_misalignment=np.array(LS_misalignment)
    Island_Piston=np.array(Island_Piston)
    
    ## call the size of the grid
    n = int(gridsize)
    
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio) # define the simualtion pupil
    lamda=proper.prop_get_wavelength(wfo) #save the wavelength value [m] into lamda
    print('beam_ratio',beam_ratio)
    if (Debug_print==True):
        print("lambda: ", lamda)

    pupil(wfo, CAL, npupil, diam, r_obstr, spiders_width, spiders_angle, pupil_file, missing_segments_number, Debug, Debug_print)

    if (Debug==True):
        fits.writeto(path + prefix+'_pupil_pre_define.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)

    proper.prop_define_entrance(wfo) #define the entrance wavefront

#wfo.wfarr *= 1./np.amax(wfo._wfarr) # max(amplitude)=1


    if (isinstance(atm_screen, (list, tuple, np.ndarray)) == True) and (atm_screen.ndim >= 2): # when the atmosphere is present
        print('atmosphere')
        atmosphere(wfo, npupil, atm_screen, Debug_print, Debug)
    
    if (isinstance(NCPA, (list, tuple, np.ndarray)) == True) and (NCPA.ndim >= 2): # when the atmosphere is present
        NCPA_application(wfo, npupil, NCPA, path, Debug_print, Debug)
    
    if (RAVC == True) or (isinstance(phase_apodizer_file, (list, tuple, np.ndarray)) == True) or (isinstance(amplitude_apodizer_file, (list, tuple, np.ndarray)) == True): # when tha apodizer is present
        apodization(wfo, r_obstr, npupil, RAVC, phase_apodizer_file, amplitude_apodizer_file, apodizer_misalignment, Debug_print, Debug)
    
    if (all(v == 0 for v in Island_Piston) == False): # when the piston is present
        island_effect_piston(wfo, npupil, Island_Piston, path, Debug_print, Debug)
    
    if (TILT.any != 0.): # when tip/tilt
        if (Debug_print==True):
            print("TILT: ", TILT)
            print("lamda: ", lamda)
        tiptilt=(np.multiply(TILT, lamda))/4 # translate the tip/tilt from lambda/D into RMS phase errors
        proper.prop_zernikes(wfo, [2,3], tiptilt) # 2-->xtilt, 3-->ytilt

    if (Debug==True):
        if CAL==1:
            fits.writeto(path + prefix+'_pupil_amplitude_CAL1.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
            fits.writeto(path + prefix+'_pupil_phase_CAL1.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        else:
            fits.writeto(path + prefix+'_pupil_amplitude_CAL0_RA'+str(int(RAVC))+'_charge'+str(charge)+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
            fits.writeto(path + prefix+'_pupil_phase_CAL0_RA'+str(int(RAVC))+'_charge'+str(charge)+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)


           
    proper.prop_propagate(wfo, f_lens, 'inizio') # propagate wavefront

    proper.prop_lens(wfo, f_lens, 'focusing lens vortex') # propagate through a lens
    proper.prop_propagate(wfo, f_lens, 'VC') # propagate wavefront

    vortex(wfo, CAL, charge, f_lens, path, Debug_print)



    if (Debug==True):
        if CAL==1:
            fits.writeto(path + prefix+'_afterVortex_CAL1.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
            fits.writeto(path + prefix+'_afterVortex_CAL1_phase.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        else:
            print('ATM: ', str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True)))
            if ((((int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))))==1):
                print('atm_screen: ', atm_screen.shape)
                print('ATM: ', int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))
            fits.writeto(path + prefix+'_afterVortex_CAL0_RA'+str(int(RAVC))+'_charge'+str(charge)+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
            fits.writeto(path + prefix+'_afterVortex_phase_CAL0_RA'+str(int(RAVC))+'_charge'+str(charge)+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_phase(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)

    proper.prop_propagate(wfo, f_lens, 'Lyot Collimetor') # propagate wavefront

    proper.prop_lens(wfo, f_lens, 'Lyot Collimetor') # propagate wavefront through  a lens
    proper.prop_propagate(wfo, f_lens, 'Lyot Stop') # propagate wavefront

    if (Debug==True):
        if CAL==1:
            fits.writeto(path + prefix+'_beforeLS_CAL1.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        else:
            fits.writeto(path + prefix+'_beforeLS_CAL0_charge'+str(charge)+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)

    lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, path, Debug_print, Debug)

    if (Debug==True):
        if CAL==1:
            fits.writeto(path + prefix+'_afterLS_CAL1.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
        else:
            fits.writeto(path + prefix+'_afterLS_CAL0_charge'+str(charge)+'_LS'+str(int(LS))+'_RA'+str(int(RAVC))+'_ATM'+str(int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))+'.fits', proper.prop_get_amplitude(wfo)[int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50),int(n/2)-int(npupil/2 + 50):int(n/2)+int(npupil/2 + 50)], overwrite=True)
 

    proper.prop_propagate(wfo, f_lens) # propagate wavefront
    proper.prop_lens(wfo, f_lens) # propagate wavefront through a lens
    proper.prop_propagate(wfo, f_lens) # propagate wavefront
    
    (wfo, sampling) = proper.prop_end(wfo, NOABS = True) # conclude the simulation --> noabs= the wavefront array will be complex
    
    
    return (wfo, sampling) # return the wavefront




