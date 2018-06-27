from __future__ import absolute_import

import numpy as np

import proper
import math

from astropy.io import fits
import os
import sys

PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))
sys.path.append(PACKAGE_PATH)

def metiscoronagraphsimulator(n,lam, pixelsize, prefix, path, diam=37., r_obstr=0.3, f_lens=658.6, pupil_file=0,  spiders_width=0.60, spiders_angle=[0., 60., 120.], charge=0, LS_parameters=[0.0, 0.0, 0.0], amplitude_apodizer_file=0,phase_apodizer_file=0,LS_amplitude_apodizer_file=0,LS_phase_apodizer_file=0,  TILT=[0.0, 0.0], atm_screen=0., missing_segments_number=0, apodizer_misalignment=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], LS_misalignment=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], Island_Piston=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], NCPA=0., NoCoro_psf=0, Offaxis_psf=0, ELT_circ=True, Vortex=False, Back=False, RAVC=False, LS=False, Debug=True, Debug_print=True, Norm_max=True, Norm_flux1=False):


    TILT=np.array(TILT)
    apodizer_misalignment=np.array(apodizer_misalignment)
    LS_misalignment=np.array(LS_misalignment)
    Island_Piston=np.array(Island_Piston)


    if (isinstance(NoCoro_psf, (list, tuple, np.ndarray)) != True):
## Non coronagraphic PSF --> simulate a non-coronagraphic psf
        (wfo_noCoro, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix, 'path':path, 'charge':0, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':0, 'amplitude_apodizer_file':0, 'TILT':[0.,0.],'LS':False,'RAVC':False, 'LS_phase_apodizer_file':0, 'LS_amplitude_apodizer_file':0,'LS_parameters':[0., 0., 0.], 'atm_screen':0, 'missing_segments_number':0, 'apodizer_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'LS_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Island_Piston':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'NCPA':0, 'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
        NoCoro_psf = (abs(wfo_noCoro))**2
        fits.writeto(path+prefix+'_psf_noCoro_nonorm.fits', NoCoro_psf[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        #if ("Norm_max" in kwargs and kwargs["Norm_max"]):
        if (Norm_max == True):
            psf_noCoro_maxnorm = NoCoro_psf/np.max(NoCoro_psf)
            fits.writeto(path+prefix+'_psf_noCoro_maxnorm.fits', psf_noCoro_maxnorm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        #if ("Norm_flux1" in kwargs and kwargs["Norm_flux1"]):
        if (Norm_flux1 == True):
            psf_noCoro_flux1norm = NoCoro_psf/sum(sum(NoCoro_psf))
            fits.writeto(path+prefix+'_psf_noCoro_flux1norm.fits', psf_noCoro_flux1norm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

    ## Coronagraphic PSF
    if (Vortex == True):
    
        if (Back == True):
            
        ## A/R --> simulate a perfect vortex by propagating a perfectly circular pupil through the vortex to the Lyot Stop, null the amplitude inside (as theory requires), then propagating back to the vortex level and save a "modified" vortex, to use in the future simulations
            (wfo_AR, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix, 'path':path, 'charge':charge, 'CAL':1, 'diam':diam, 'spiders_width':0, 'spiders_angle':[0., 0., 0.], 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':0., 'pupil_file':0, 'phase_apodizer_file':0, 'amplitude_apodizer_file':0, 'TILT':[0.,0.],'LS':False,'RAVC':False, 'LS_phase_apodizer_file':0, 'LS_amplitude_apodizer_file':0,'LS_parameters':[0., 0., 0.], 'atm_screen':0, 'missing_segments_number':0, 'apodizer_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'LS_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Island_Piston':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'NCPA':0, 'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
        
        if (isinstance(Offaxis_psf, (list, tuple, np.ndarray)) != True):
            ## No Vortex Si Mask --> simulate a non-coronagraphic-apodized psf: apodizer and Lyot Stop are present, but not the vortex --> as an off-axis psf
            (wfo_offaxis, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix, 'path':path, 'charge':0, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':phase_apodizer_file, 'amplitude_apodizer_file':amplitude_apodizer_file, 'TILT':[0.,0.],'LS':LS,'RAVC':RAVC, 'LS_phase_apodizer_file':LS_phase_apodizer_file, 'LS_amplitude_apodizer_file':LS_amplitude_apodizer_file,'LS_parameters':LS_parameters, 'atm_screen':0, 'missing_segments_number':missing_segments_number, 'apodizer_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'LS_misalignment':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Island_Piston':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'NCPA':0, 'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
            Offaxis_psf = (abs(wfo_offaxis))**2
            fits.writeto(path+prefix+'_psf_offaxis_nonorm.fits', Offaxis_psf[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
            if (Norm_max == True):
                #if ("Norm_max" in kwargs and kwargs["Norm_max"]):
                psf_noVortex_Mask_maxnorm = Offaxis_psf/np.max(NoCoro_psf)
                fits.writeto(path+prefix+'_psf_offaxis_maxnorm.fits', psf_noVortex_Mask_maxnorm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
            if (Norm_flux1 == True):
                #if ("Norm_flux1" in kwargs and kwargs["Norm_flux1"]):
                psf_noVortex_Mask_flux1norm = Offaxis_psf/sum(sum(NoCoro_psf))
                fits.writeto(path+prefix+'_psf_offaxis_flux1norm.fits', psf_noVortex_Mask_flux1norm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

    atm_screen=np.array(atm_screen)
    NCPA=np.array(NCPA)

    if (atm_screen.ndim == 3) or (TILT.ndim == 2) or (LS_misalignment.ndim == 2) or (apodizer_misalignment.ndim == 2) or (Island_Piston.ndim == 2) or (NCPA.ndim == 3):
        print('Cube')
        
        if (atm_screen.ndim == 3):
            length_cube = atm_screen.shape[0]
        if (TILT.ndim == 2):
            length_cube = TILT.shape[0]
        if (LS_misalignment.ndim == 2):
            length_cube = LS_misalignment.shape[0]
        if (apodizer_misalignment.ndim == 2):
            length_cube = apodizer_misalignment.shape[0]
        if (Island_Piston.ndim == 2):
            length_cube = Island_Piston.shape[0]
        if (NCPA.ndim == 3):
            length_cube = NCPA.shape[0]
        
        psf_Coro = np.zeros((length_cube,n,n))

        
        for iter in range(0, length_cube):
            print('iter: ', iter)

            if ((isinstance(atm_screen, (list, tuple, np.ndarray)) == True)):
                if (atm_screen.ndim == 3):
                    atm_screen_iter = atm_screen[iter,:,:]
                else:
                    atm_screen_iter = atm_screen
            if (TILT.ndim == 2):
                TILT_iter = TILT[iter,:]
                print('TILT: ', TILT_iter)
            else:
                TILT_iter = TILT
            if (LS_misalignment.ndim == 2):
                LS_misalignment_iter =  LS_misalignment[iter,:]
            else:
                LS_misalignment_iter =  LS_misalignment
            if (apodizer_misalignment.ndim == 2):
                apodizer_misalignment_iter = apodizer_misalignment[iter,:]
            else:
                apodizer_misalignment_iter = apodizer_misalignment
            if (Island_Piston.ndim == 2):
                Island_Piston_iter = Island_Piston[iter,:]
            else:
                Island_Piston_iter = Island_Piston
            if (isinstance(NCPA, (list, tuple, np.ndarray)) == True):
                if (NCPA.ndim == 3):
                    NCPA_iter = NCPA[iter,:,:]
                else:
                    NCPA_iter = NCPA

            if (Vortex == True):
                ## Si Vortex Si Mask --> simulate the coronagraphic-apodized psf
                (wfo_Coro, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix,'path':path,  'charge':charge, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':phase_apodizer_file, 'amplitude_apodizer_file':amplitude_apodizer_file, 'TILT':TILT_iter,'LS':LS,'RAVC':RAVC, 'LS_phase_apodizer_file':LS_phase_apodizer_file, 'LS_amplitude_apodizer_file':LS_amplitude_apodizer_file,'LS_parameters':LS_parameters,  'atm_screen':atm_screen_iter, 'missing_segments_number':missing_segments_number, 'apodizer_misalignment':apodizer_misalignment_iter, 'LS_misalignment':LS_misalignment_iter, 'Island_Piston':Island_Piston_iter, 'NCPA':NCPA_iter,'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
                psf_Coro[iter,:,:] = (abs(wfo_Coro))**2
            else:
        
            ## Apodizer --> simulate a coronagraphic psf with an apodizer (no Vortex)
                (wfo_apodizer, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix,'path':path,  'charge':0, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':phase_apodizer_file, 'amplitude_apodizer_file':amplitude_apodizer_file, 'TILT':TILT_iter,'LS':LS,'RAVC':False, 'LS_phase_apodizer_file':LS_phase_apodizer_file, 'LS_amplitude_apodizer_file':LS_amplitude_apodizer_file,'LS_parameters':LS_parameters, 'atm_screen':atm_screen_iter, 'missing_segments_number':missing_segments_number, 'apodizer_misalignment':apodizer_misalignment_iter, 'LS_misalignment':LS_misalignment_iter, 'Island_Piston':Island_Piston_iter,'NCPA':NCPA_iter, 'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
                psf_Coro[iter,:,:]  = (abs(wfo_apodizer))**2
        
        fits.writeto(path+prefix+'_psf_cube_Coro_nonorm.fits', psf_Coro[:,int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        if (Norm_max == True):
            if (Vortex == True):
                psf_Coro_maxnorm = psf_Coro/np.max(Offaxis_psf)
            else:
                psf_Coro_maxnorm = psf_Coro/np.max(NoCoro_psf)
        fits.writeto(path+prefix+'_psf_cube_Coro_maxnorm.fits', psf_Coro_maxnorm[:,int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        if (Norm_flux1 == True):
            psf_Coro_flux1norm = psf_Coro/sum(sum(NoCoro_psf))
            fits.writeto(path+prefix+'_psf_cube_Coro_flux1norm.fits', psf_Coro_flux1norm[:,int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

        
    else:
        print('No Cube')
        if (Vortex == True):
            print('VC') 
            print('ATM: ', int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))
            print('atm screen: ', atm_screen.shape)
        ## Si Vortex Si Mask --> simulate the coronagraphic-apodized psf
            (wfo_Coro, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix,'path':path,  'charge':charge, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':phase_apodizer_file, 'amplitude_apodizer_file':amplitude_apodizer_file, 'TILT':TILT,'LS':LS,'RAVC':RAVC, 'LS_phase_apodizer_file':LS_phase_apodizer_file, 'LS_amplitude_apodizer_file':LS_amplitude_apodizer_file,'LS_parameters':LS_parameters,  'atm_screen':atm_screen, 'missing_segments_number':missing_segments_number, 'apodizer_misalignment':apodizer_misalignment, 'LS_misalignment':LS_misalignment, 'Island_Piston':Island_Piston, 'NCPA':NCPA,'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
            psf_Coro = (abs(wfo_Coro))**2
        
        else:
            print('Apodizer')  
            print('ATM: ', int(isinstance(atm_screen, (list, tuple, np.ndarray)) == True))
            print('atm screen: ', atm_screen.shape)
            ## Apodizer --> simulate a coronagraphic psf with an apodizer (no Vortex)
            (wfo_apodizer, sampling) = proper.prop_run('telescope', lam, n, PASSVALUE={'prefix':prefix,'path':path,  'charge':0, 'CAL':0, 'diam':diam, 'spiders_width':spiders_width, 'spiders_angle':spiders_angle, 'beam_ratio': beam_ratio, 'f_lens':f_lens, 'npupil':npupil, 'r_obstr':r_obstr, 'pupil_file':pupil_file, 'phase_apodizer_file':phase_apodizer_file, 'amplitude_apodizer_file':amplitude_apodizer_file, 'TILT':TILT,'LS':LS,'RAVC':False, 'LS_phase_apodizer_file':LS_phase_apodizer_file, 'LS_amplitude_apodizer_file':LS_amplitude_apodizer_file,'LS_parameters':LS_parameters, 'atm_screen':atm_screen, 'missing_segments_number':missing_segments_number, 'apodizer_misalignment':apodizer_misalignment, 'LS_misalignment':LS_misalignment, 'Island_Piston':Island_Piston,'NCPA':NCPA, 'Debug_print':Debug_print,'Debug':Debug}, QUIET=True)
            psf_Coro = (abs(wfo_apodizer))**2

        fits.writeto(path+prefix+'_psf_Coro_nonorm.fits', psf_Coro[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        if (Norm_max == True):
            if (Vortex == True):
                psf_Coro_maxnorm = psf_Coro/np.max(Offaxis_psf)
            else:
                psf_Coro_maxnorm = psf_Coro/np.max(NoCoro_psf)
            fits.writeto(path+prefix+'_psf_Coro_maxnorm.fits', psf_Coro_maxnorm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
        if (Norm_flux1 == True):
            psf_Coro_flux1norm = psf_Coro/sum(sum(NoCoro_psf))
            fits.writeto(path+prefix+'_psf_Coro_flux1norm.fits', psf_Coro_flux1norm[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)



    return

if __name__ == '__main__':
    metiscoronagraphsimulator()

    
