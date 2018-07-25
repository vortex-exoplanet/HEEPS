#!/usr/bin/env python3

from heeps import apodization, vortex,lyotstop # loads all HEEPS scripts required for simuation


def coronagraphs(wfo, r_obstr,npupil, phase_apodizer_file,amplitude_apodizer_file,apodizer_misalignment,charge,f_lens,diam,LS_amplitude_apodizer_file,LS_misalignment,LS,LS_parameters,spiders_angle, LS_phase_apodizer_file, Debug_print,pixelsize, Debug, coronagraph_type='None'):
    
    if coronagraph_type == 'RAVC':
        phase_apodizer_file = 0
        RAVC = True
        apodization(wfo, r_obstr, npupil, RAVC=True, phase_apodizer_file=phase_apodizer_file, amplitude_apodizer_file=amplitude_apodizer_file, apodizer_misalignment=apodizer_misalignment, Debug_print=Debug_print)
        vortex(wfo, charge, f_lens,diam, pixelsize, Debug_print = Debug_print)
        lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug)

    elif coronagraph_type == "VC":
        phase_apodizer_file = 0
        RAVC = False
        vortex(wfo, charge, f_lens,diam, pixelsize, Debug_print = Debug_print)
        lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug)

    elif coronagraph_type == 'APP':
        RAVC = False
        charge = 0
        apodization(wfo, r_obstr, npupil, RAVC=False, phase_apodizer_file=phase_apodizer_file, amplitude_apodizer_file=amplitude_apodizer_file, apodizer_misalignment=apodizer_misalignment, Debug_print=Debug_print)
        vortex(wfo, charge, f_lens,diam, pixelsize, Debug_print = Debug_print)
        lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug)

    elif coronagraph_type == 'OFFAXIS':
        print('No Coronagraph')    
        RAVC = False
        phase_apodizer_file = 0
        charge = 0
#        apodization(wfo, r_obstr, npupil, RAVC=False, phase_apodizer_file=phase_apodizer_file, amplitude_apodizer_file=amplitude_apodizer_file, apodizer_misalignment=apodizer_misalignment, Debug_print=Debug_print)
#        vortex(wfo, charge, f_lens,diam, pixelsize, Debug_print = Debug_print)
        lyotstop(wfo, diam, r_obstr, npupil, RAVC, LS, LS_parameters, spiders_angle, LS_phase_apodizer_file, LS_amplitude_apodizer_file, LS_misalignment, Debug_print, Debug)
    else:
        print('ELT PSF')    
    return wfo


 