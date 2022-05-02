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

    if mode in ['CVC', 'RAVC', 'CLC']:

        # case 1: load mask from data
        if ls_mask is not None:
            if verbose is True:
                print("   apply lyot stop data from 'ls_mask'")
            ls_mask = resize_img(ls_mask, npupil)

        # case 2: load mask from file
        elif os.path.isfile(f_lyot_stop):
            if verbose is True:
                print("   apply lyot stop from '%s'"%os.path.basename(f_lyot_stop))
            ls_mask = resize_img(fits.getdata(f_lyot_stop), npupil)

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
                print('   apply Lyot stop: circ_ext=%s, circ_int=%s'
                    %(circ_ext, circ_int)
                    + ', ls_dRext=%.4f, ls_dRint=%.4f, ls_dRspi=%.4f'
                    %(ls_dRext, ls_dRint, ls_dRspi))
            
        # optional: can return updated conf, instead of applying Lyot stop
        if not apply_ls:
            return ls_mask
        else:
            proper.prop_multiply(wf, pad_img(ls_mask, ngrid))

    return wf