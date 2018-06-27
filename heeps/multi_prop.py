from __future__ import absolute_import

import numpy as np
import proper
import math
from astropy.io import fits
#from .telescope import telescope as telescope
import os
import sys

PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))
sys.path.append(PACKAGE_PATH)

def multi_prop(n,lam, prefix, beam_ratio, RAVC=True, LS=True):
    
    diam = 37.
    print("lam: ", lam)
    print("beam_ratio: ", beam_ratio)
    npupil = math.ceil(n*beam_ratio) # compute the pupil size --> has to be ODD (proper puts the center in the up right pixel next to the grid center)
    if npupil % 2 == 0:
        npupil = npupil +1


    
        ## Non coronagraphic PSF --> simulate a non-coronagraphic psf
    (wfo_noCoro, sampling) = proper.prop_run_multi('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':0, 'diam':diam,'beam_ratio':beam_ratio, 'RAVC':False, 'LS':False}, QUIET=True)
#(wfo_noCoro, sampling) = proper.prop_run('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':0, 'diam':37., 'RAVC':False, 'LS':False}, QUIET=True)
    NoCoro_psf = (abs(wfo_noCoro))**2
    fits.writeto(prefix+'psf_noCoro_nonorm.fits', NoCoro_psf, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_noCoro_maxnorm = NoCoro_psf/np.max(NoCoro_psf)
    fits.writeto(prefix+'_psf_noCoro_maxnorm.fits', psf_noCoro_maxnorm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_noCoro_flux1norm = NoCoro_psf/sum(sum(NoCoro_psf))
    fits.writeto(prefix+'_psf_noCoro_flux1norm.fits', psf_noCoro_flux1norm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

    ## No Vortex Si Mask --> simulate a non-coronagraphic-apodized psf: apodizer and Lyot Stop are present, but not the vortex --> as an off-axis psf
    (wfo_offaxis, sampling) = proper.prop_run_multi('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':0, 'diam':diam,'beam_ratio':beam_ratio, 'RAVC':RAVC, 'LS':LS}, QUIET=True)
#(wfo_offaxis, sampling) = proper.prop_run('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':0, 'diam':37., 'RAVC':RAVC, 'LS':LS}, QUIET=True)
    Offaxis_psf = (abs(wfo_offaxis))**2
    fits.writeto(prefix+'psf_offaxis_nonorm.fits', Offaxis_psf, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_noVortex_Mask_maxnorm = Offaxis_psf/np.max(NoCoro_psf)
    fits.writeto(prefix+'psf_offaxis_maxnorm.fits', psf_noVortex_Mask_maxnorm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_noVortex_Mask_flux1norm = Offaxis_psf/sum(sum(NoCoro_psf))
    fits.writeto(prefix+'psf_offaxis_flux1norm.fits', psf_noVortex_Mask_flux1norm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

    (wfo_Coro, sampling) = proper.prop_run_multi('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':2, 'diam':diam,'beam_ratio':beam_ratio, 'RAVC':RAVC, 'LS':LS}, QUIET=True)
#(wfo_Coro, sampling) = proper.prop_run('simple_telescope', lam, n, PASSVALUE={'prefix':prefix, 'charge':2, 'diam':37., 'RAVC':RAVC, 'LS':LS}, QUIET=True)
    psf_Coro = (abs(wfo_Coro))**2
    fits.writeto(prefix+'_psf_Coro_nonorm.fits', psf_Coro, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_Coro_maxnorm = psf_Coro/np.max(Offaxis_psf)
    fits.writeto(prefix+'_psf_Coro_maxnorm.fits', psf_Coro_maxnorm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)
    psf_Coro_flux1norm = psf_Coro/sum(sum(NoCoro_psf))
    fits.writeto(prefix+'_psf_Coro_flux1norm.fits', psf_Coro_flux1norm, overwrite=True)#[int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio),int(n/2)-math.ceil(50./beam_ratio):int(n/2)+math.ceil(50./beam_ratio)], overwrite=True)

    

    return

