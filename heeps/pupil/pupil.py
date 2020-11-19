from .create_pupil import create_pupil
from heeps.util.img_processing import resize_img, pad_img
from heeps.util.save2fits import save2fits
import proper
import numpy as np
import os.path
from astropy.io import fits 


def pupil(file_pupil='', lam=3.8e-6, ngrid=1024, npupil=285, pupil_img_size=40, diam_ext=37, 
        diam_int=11, spi_width=0.5, spi_angles=[0,60,120], npetals=6, seg_width=0, 
        seg_gap=0, seg_rms=0, seg_ny=[10,13,16,19,22,23,24,25,26,27,28,29,
        30,31,30,31,30,31,30,31,30,31,30,29,28,27,26,25,24,23,22,19,16,13,10],
        seg_missing=[], select_petal=None, savefits=False, verbose=False, **conf):
    
    ''' Create a wavefront object at the entrance pupil plane. 
    The pupil is either loaded from a fits file, or created using 
    pupil parameters.
    Can also select only one petal and mask the others.

    Args:
        dir_output (str):
            path to saved pupil file
        band (str):
            spectral band (e.g. 'L', 'M', 'N1', 'N2')
        mode (str):
            HCI mode: RAVC, CVC, APP, CLC
        file_pupil: str
            path to a pupil fits file
        lam: float
            wavelength in m
        ngrid: int
            number of pixels of the wavefront array
        npupil: int
            number of pixels of the pupil
        pupil_img_size: float
            pupil image (for PROPER) in m
        diam_ext: float
            outer circular aperture in m
        diam_int: float
            central obscuration in m
        spi_width: float
            spider width in m
        spi_angles: list of float
            spider angles in deg
        seg_width: float
            segment width in m
        seg_gap: float
            gap between segments in m
        seg_rms: float
            rms of the reflectivity of all segments
        seg_ny: list of int
            number of hexagonal segments per column (from left to right)
        seg_missing: list of tupples
            coordinates of missing segments
        npetals: int
            number of petals
        select_petal: int
            selected petal in range npetals, default to None
    
    '''

    # initialize wavefront using PROPER
    beam_ratio = npupil/ngrid*(diam_ext/pupil_img_size)
    wf = proper.prop_begin(pupil_img_size, lam, ngrid, beam_ratio)

    # load pupil file
    if os.path.isfile(file_pupil):
        if verbose is True:
            print("Load pupil from '%s'"%os.path.basename(file_pupil))
        #conf.update()
        pup = fits.getdata(file_pupil).astype(np.float32)
        # resize to npupil
        pup = resize_img(pup, npupil)
    # if no pupil file, create a pupil
    else:
        if verbose is True:
            print("Create pupil: spi_width=%s m, seg_width=%s m, seg_gap=%s m, seg_rms=%s"\
                %(spi_width, seg_width, seg_gap, seg_rms))
        conf.update(npupil=npupil,
                    pupil_img_size=pupil_img_size, 
                    diam_ext=diam_ext, 
                    diam_int=diam_int, 
                    spi_width=spi_width, 
                    spi_angles=spi_angles, 
                    seg_width=seg_width, 
                    seg_gap=seg_gap, 
                    seg_ny=seg_ny, 
                    seg_missing=seg_missing,
                    seg_rms=seg_rms)
        pup = create_pupil(**conf)

    # normalize the entrance pupil intensity (total flux = 1)
    I_pup = pup**2
    pup = np.sqrt(I_pup/np.sum(I_pup))
    # pad with zeros and add to wavefront
    proper.prop_multiply(wf, pad_img(pup, ngrid))
        
    # select one petal (optional)
    if select_petal in range(npetals) and npetals > 1:
        if verbose is True:
            print("   select_petal=%s"%select_petal)
        # petal start and end angles
        pet_angle = 2*np.pi/npetals
        pet_start = pet_angle/2 + (select_petal - 1)*pet_angle
        pet_end = pet_start + pet_angle
        # petal angles must be 0-2pi
        ti %= (2*np.pi)
        pet_start %= (2*np.pi)
        pet_end %= (2*np.pi)
        # check if petal crosses 0 angle
        if pet_end - pet_start < 0:
            pet_start = (pet_start + np.pi)%(2*np.pi) - np.pi
            pet_end = (pet_end + np.pi)%(2*np.pi) - np.pi
            ti = (ti + np.pi)%(2*np.pi) - np.pi
        # create petal and add to grid
        petal = ((ti>=pet_start) * (ti<=pet_end)).astype(int)
        petal = resize_img(petal, ngrid)
        proper.prop_multiply(wf, petal)

    # save pupil as fits file
    if savefits == True:
        save2fits(pup, 'pupil', **conf)

    if verbose is True:
        print('   diam=%s m, resize to %s pix, zero-pad to %s pix\n'\
            %(round(diam_ext, 2), npupil, ngrid))
    
    return wf