import heeps.util.img_processing as impro
import proper
import numpy as np

def static_ncpa(wf, ncpa_screen, **conf):
    
    # get useful parameters
    gridsize = conf['gridsize']
    npupil = conf['npupil']
    
    # resize the phase screen, and pad with zeros to match PROPER gridsize
    ncpa_screen = impro.resize_img(ncpa_screen, npupil)
    ncpa_screen = impro.pad_img(ncpa_screen, gridsize)
    # add the phase screen
    proper.prop_add_phase(wf, ncpa_screen)
    
    return wf
