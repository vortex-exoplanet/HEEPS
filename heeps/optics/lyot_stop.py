from heeps.util.img_processing import resize_img, pad_img
import proper
from astropy.io import fits
import os.path

def lyot_stop(wf, ls_mask=None, f_lyot_stop='', ngrid=1024, npupil=285, 
        mode='RAVC', ravc_r=0.6, ls_dRext=0.03, ls_dRint=0.05, ls_dRspi=0.04, 
        spi_width=0.5, spi_angles=[0,60,120], diam_ext=37, diam_int=11, 
        diam_nominal=37, ls_misalign=None, verbose=False, **conf):

    """ Add a Lyot stop, in the case of a focal plane mask. """
    
    if mode in ['CVC', 'RAVC', 'CLC']:

        # case 1: load mask from data
        if ls_mask is not None:
            if verbose is True:
                print("   apply lyot stop data from 'ls_mask'")
            ls_mask = resize_img(ls_mask, npupil)
            proper.prop_multiply(wf, pad_img(ls_mask, ngrid))

        # case 2: load mask from file
        elif os.path.isfile(f_lyot_stop):
            if verbose is True:
                print("   apply lyot stop from '%s'"%os.path.basename(f_lyot_stop))
            ls_mask = resize_img(fits.getdata(f_lyot_stop), npupil)
            proper.prop_multiply(wf, pad_img(ls_mask, ngrid))

        # case 3: create a lyot stop mask
        else:
            # scale nominal values to pupil external diameter
            scaling = diam_nominal/diam_ext
            # LS parameters
            r_obstr = ravc_r if mode in ['RAVC'] else diam_int/diam_ext
            ls_int = r_obstr + ls_dRint*scaling
            ls_ext = 1 - ls_dRext*scaling
            ls_spi = spi_width/diam_ext + ls_dRspi*scaling
            # LS misalignments
            ls_misalign = [0,0,0,0,0,0] if ls_misalign is None else list(ls_misalign)
            dx_amp, dy_amp = ls_misalign[0:2]
            # create Lyot stop
            proper.prop_circular_aperture(wf, ls_ext, dx_amp, dy_amp, NORM=True)
            if diam_int > 0:
                proper.prop_circular_obscuration(wf, ls_int, dx_amp, dy_amp, NORM=True)
            if spi_width > 0:
                for angle in spi_angles:
                    proper.prop_rectangular_obscuration(wf, 2*ls_spi, 2, 
                            dx_amp, dy_amp, ROTATION=angle, NORM=True)
            if verbose is True:
                print('   apply Lyot stop: ls_int=%s, ls_ext=%s, ls_spi=%s'
                    %(round(ls_int, 4), round(ls_ext, 4), round(ls_spi, 4)))

    return wf