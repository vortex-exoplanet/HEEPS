import heeps.util.img_processing as impro
import proper
import numpy as np
from heeps.coronagraphs import circular_apodization


def apodization(wf, conf, margin=0, RAVC=False, RAVC_amp_file=None, RAVC_phase_file=None):
    
    if RAVC is False:
        return wf, None, None
    else:
        # load parameters
        gridsize = conf['gridsize']
        npupil = conf['npupil']
        R_obstr = conf['R_obstr']
        get_amp = conf['get_amp']
        get_phase = conf['get_phase']
        
        # get apodizer misalignments
        if conf['RAVC_misalign'] is not None:
            RAVC_misalign = list(conf['RAVC_misalign'])
        else:
            RAVC_misalign = [0,0,0,0,0,0]
        dx_amp, dy_amp, dz_amp = RAVC_misalign[0:3]
        dx_phase, dy_phase, dz_phase = RAVC_misalign[3:6]
        
        if None in [RAVC_amp_file, RAVC_phase_file]:
            # define the apodizer transmission and radius [Mawet2013]
            T_RAVC = 1 - (R_obstr**2 + R_obstr*np.sqrt(R_obstr**2 + 8))/4
            R_RAVC = R_obstr/np.sqrt(1 - T_RAVC)
            # define the apodizer as a ring (with % misalignments)
            apodizer = circular_apodization(wf, R_RAVC, 1., T_RAVC, \
                    xc=RAVC_misalign[0], yc=RAVC_misalign[1], NORM=True)
            apodizer = proper.prop_shift_center(apodizer)
        else:
            # get amplitude and phase files
            RAVC_amp_file = os.path.join(conf['input_dir'], conf['RAVC_amp_file'])
            RAVC_phase_file = os.path.join(conf['input_dir'], conf['RAVC_phase_file'])
            # get amplitude and phase data
            RAVC_amp = fits.getdata(RAVC_amp_file) if os.path.isfile(RAVC_amp_file) \
                    else np.ones((npupil, npupil))
            RAVC_phase = fits.getdata(RAVC_phase_file) if os.path.isfile(RAVC_phase_file) \
                    else np.zeros((npupil, npupil))
            # resize to npupil
            RAVC_amp = impro.resize_img(RAVC_amp, npupil)
            RAVC_phase = impro.resize_img(RAVC_phase, npupil)
            # pad with zeros to match PROPER gridsize
            RAVC_amp = impro.pad_img(RAVC_amp, gridsize)
            RAVC_phase = impro.pad_img(RAVC_phase, gridsize)
            # apodizer
            apodizer = RAVC_amp*np.exp(1j*RAVC_phase)
        
        # multiply the loaded apodizer
        proper.prop_multiply(wf, apodizer)
        
        # get the apodizer amplitude and phase for output
        if conf['full_output'] is True:
            apo_amp = impro.crop_img(proper.prop_get_amplitude(wf), npupil, margin)
            apo_phase = impro.crop_img(proper.prop_get_phase(wf), npupil, margin)
            return wf, apo_amp, apo_phase
        else:
            return wf, None, None