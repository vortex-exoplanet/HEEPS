from heeps.util.img_processing import resize_img, pad_img
from heeps.pupil.create_stop import create_stop
import proper
from astropy.io import fits
import os.path

def lyot_stop(wf, ls_mask=None, f_lyot_stop='', ngrid=1024, npupil=285, 
        mode='RAVC', ravc_r=0.6, ls_dRext=0.03, ls_dRint=0.05, ls_dRspi=0.04, 
        diam_ext=37, diam_int=11, ls_ext_circ=True, ls_int_circ=True, 
        ls_misalign=None, apply_ls=True, verbose=False, **conf):

    """ Add a Lyot stop, in the case of a focal plane mask. """

    if mode in ['CVC', 'RAVC', 'CLC', 'IMG', 'LMS']:

        # case 1: mask already preloaded
        if ls_mask is not None:
            if verbose is True:
                print("   apply preloaded lyot stop")

        # case 2: load mask from file
        elif os.path.isfile(f_lyot_stop):
            if verbose is True:
                print("   apply lyot stop from '%s'"%os.path.basename(f_lyot_stop))
            ls_mask = resize_img(fits.getdata(f_lyot_stop), npupil)
            ls_mask = pad_img(ls_mask, ngrid)

        # case 3: create a lyot stop mask
        else:
            # load params
            d_ext = diam_ext
            d_int = ravc_r*diam_ext if 'RAVC' in mode else diam_int
            circ_ext = ls_ext_circ
            circ_int = True         if 'RAVC' in mode else ls_int_circ
            misalign_x , misalign_y = [0, 0] if ls_misalign is None \
                                            else list(ls_misalign)[0:2]
            # create Lyot stop
            ls_mask = create_stop(d_ext, d_int, ls_dRext, ls_dRint, ls_dRspi,
                npupil=npupil, misalign_x=misalign_x, misalign_y=misalign_y,
                circ_ext=circ_ext, circ_int=circ_int, **conf)
            if verbose is True:
                print('   apply Lyot stop: circ_ext/int=%s'%[circ_ext, circ_int]
                    + ', ls_dRext/int/spi=%s'%[ls_dRext, ls_dRint, ls_dRspi]
                    + ', ls_misalign=%s'%ls_misalign)
            # zero-pad
            ls_mask = pad_img(ls_mask, ngrid)

        # optional: can return lyot stop array, instead of applying to wavefront
        if not apply_ls:
            return ls_mask
        else:
            proper.prop_multiply(wf, ls_mask)

    return wf