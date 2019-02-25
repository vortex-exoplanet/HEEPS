import heeps.util.img_processing as impro
import proper
import numpy as np
import astropy.units as u
from astropy.io import fits
import os.path

def island_effect_piston(wf, petal_piston=[0,0,0,0,0,0], npetals=6, **conf):
    
    # get useful parameters
    lam = conf['lam']
    gridsize = conf['gridsize']
    npupil = conf['npupil']
    
    # path to petal piston files
    filename = os.path.join(conf['input_dir'], '1024_pixelsize5mas_petal%s_243px.fits')
    
    # multiply all pistons by their respective petal, and sum them up
    piston_screen = np.sum(np.float32([petal_piston[x] \
            *fits.getdata(filename%(x+1)) for x in range(npetals)]), 0)
    
    # resize and pad with zeros to match PROPER gridsize
    piston_screen = impro.resize_img(piston_screen, npupil)
    piston_screen = impro.pad_img(piston_screen, gridsize)
    
    # wavenumber (spatial angular frequency) in rad / µm
    k = 2*np.pi/(lam*u.m).to('um').value
    # multiply the wavefront by the complex phase screen
    proper.prop_multiply(wf, np.exp(1j*k*piston_screen))
    
    return wf
