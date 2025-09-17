from heeps.util.img_processing import resize_img, pad_img
from heeps.util.compressed_file import read_compressed_fits
from pathlib import Path
from heeps.pupil.create_stop import create_stop
import proper
from astropy.io import fits
import os.path
import numpy as np

def lyot_stop(wf, ls_mask=None, f_lyot_stop='', f_lyot_stop_phase='',
              ngrid=1024, npupil=285, 
        mode='RAVC', ravc_r=0.6, ls_dRext=0.03, ls_dRint=0.05, ls_dRspi=0.04, 
        diam_ext=37, diam_int=11, ls_ext_circ=True, ls_int_circ=True, 
        ls_misalign=None, apply_ls=True, verbose=False, **conf):

    """
    Apply a Lyot stop, in the case of a focal plane mask. 

    Parameters
    ----------
    wf : obj
        Wavefront object, as returned by `proper.prop_init`.
    ls_mask : array, optional
        Preloaded Lyot stop mask (overrides other options).
    f_lyot_stop : str, optional
        File from which to load Lyot stop mask.
    f_lyot_stop_phase : str, optional
        File from which to load Lyot stop phase.
    ngrid : int, optional
        Size of wavefront grid.
    npupil : int, optional
        Size of pupil.
    mode : str, optional
        Mode of operation. Options are 'RAVC', 'CVC', 'CLC', 'IMG', 'LMS'.
    ravc_r : float, optional
        Radius of apodizer in RAVC mode.
    ls_dRext : float, optional
        LS Rext undersize (fraction of diam_ext).
    ls_dRint : float, optional
        LS Rint oversize (fraction of diam_ext).
    ls_dRspi : float, optional
        LS Rspi oversize (fraction of diam_ext).
    diam_ext : float, optional
        effective outer circular aperture in m
    diam_int : float, optional
        effective central obscuration in m
    ls_ext_circ : bool, optional
        If True, make outer ring of Lyot stop circular.
    ls_int_circ : bool, optional
        If True, make inner ring of Lyot stop circular.
    ls_misalign : list, optional
        Misalignment of Lyot stop (dx, dy).
    apply_ls : bool, optional
        If False, return Lyot stop mask instead of applying to wavefront.
    verbose : bool, optional
        If True, print out what's happening.

    Returns
    -------
    wf : obj
        Modified wavefront object.
    """

    if mode in ['CVC', 'RAVC', 'CLC', 'IMG', 'LMS']:

        # case 1: mask already preloaded
        if ls_mask is not None:
            if verbose is True:
                print("   apply preloaded lyot stop")

        # case 2: load mask from file
        elif os.path.isfile(f_lyot_stop):
            if verbose is True:
                print("   apply lyot stop from '%s'"%os.path.basename(f_lyot_stop))
            file_path = Path(f_lyot_stop)
            extension = "".join(file_path.suffixes)
            if extension == '.fits':
                lyot_array = fits.getdata(f_lyot_stop)
            elif extension == '.fits.tar.gz':
                lyot_array = read_compressed_fits(f_lyot_stop)
            else:
                raise ValueError('Lyot stop file must be in FITS format (.fits or .fits.tar.gz)')          
            ls_mask = resize_img(lyot_array, npupil)
            ls_mask = pad_img(ls_mask, ngrid)
        # case 3: create a lyot stop mask
        else:
            # load params
            circ_ext = ls_ext_circ
            circ_int = ls_int_circ
            if 'RAVC' in mode:
                diam_int = ravc_r*diam_ext
                circ_int = True
            dx , dy = [0, 0] if ls_misalign is None else list(ls_misalign)[0:2]
            # create Lyot stop
            ls_mask = create_stop(diam_ext=diam_ext, diam_int=diam_int,
                dRext=ls_dRext, dRint=ls_dRint, dRspi=ls_dRspi,
                circ_ext=circ_ext, circ_int=circ_int, dx=dx, dy=dy,
                npupil=npupil, **conf)
            if verbose is True:
                print('   apply Lyot stop: circ_ext/int=%s'%[circ_ext, circ_int]
                    + ', ls_dRext/int/spi=%s'%[ls_dRext, ls_dRint, ls_dRspi]
                    + ', ls_misalign=%s'%ls_misalign)
            # zero-pad
            ls_mask = pad_img(ls_mask, ngrid)
        
        if os.path.isfile(f_lyot_stop_phase):
            ls_phase = resize_img(fits.getdata(f_lyot_stop_phase), npupil)
            ls_phase = pad_img(ls_phase, ngrid)
        else:
            ls_phase = None

        # optional: can return lyot stop array, instead of applying to wavefront

        if ls_phase is not None:
            if not apply_ls:
                return ls_mask * np.exp(1j* 2*np.pi/wf.lamda * ls_phase)
            else:
                proper.prop_multiply(wf, ls_mask * np.exp(1j* 2*np.pi/wf.lamda * ls_phase))
        else:
            if not apply_ls:
                return ls_mask
            else:
                proper.prop_multiply(wf, ls_mask)

    return wf