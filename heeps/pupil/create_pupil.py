from .create_hexagon import create_hexagon
from heeps.util.img_processing import resize_img, pad_img
import proper
import numpy as np
import os.path
from astropy.io import fits 


def create_pupil(nhr=2**10, npupil=285, pupil_img_size=40, diam_ext=37, 
        diam_int=11, spi_width=0.5, spi_angles=[0,60,120], seg_width=0, 
        seg_gap=0, seg_rms=0, seg_ny=[10,13,16,19,22,23,24,25,26,27,28,29,
        30,31,30,31,30,31,30,31,30,31,30,29,28,27,26,25,24,23,22,19,16,13,10], 
        seg_missing=[], seed=123456, **conf):
    
    ''' Create a pupil.
    
    Args:
        nhr: int
            high resolution grid
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
    
    '''

    # create a high res pupil with PROPER of even size (nhr)
    nhr_size = pupil_img_size*nhr/(nhr-1)
    wf_tmp = proper.prop_begin(nhr_size, 1, nhr, diam_ext/nhr_size)
    if diam_ext > 0:
        proper.prop_circular_aperture(wf_tmp, 1, NORM=True)
    if diam_int > 0:
        proper.prop_circular_obscuration(wf_tmp, diam_int/diam_ext, NORM=True)
    if spi_width > 0:
        for angle in spi_angles:
            proper.prop_rectangular_obscuration(wf_tmp, spi_width/nhr_size, 2, \
                ROTATION=angle, NORM=True)
    pup = proper.prop_get_amplitude(wf_tmp)
    # crop the pupil to odd size (nhr-1), and resize to npupil
    pup = pup[1:,1:]
    pup = resize_img(pup, npupil)
    # add segments
    if seg_width > 0:
        segments = np.zeros((nhr, nhr))
        # sampling in meters/pixel
        sampling = pupil_img_size/nhr
        # dist between center of two segments, side by side
        seg_d = seg_width*np.cos(np.pi/6) + seg_gap
        # segment radius
        seg_r = seg_width/2
        # segment radial distance wrt x and y axis
        seg_ny = np.array(seg_ny)
        seg_nx = len(seg_ny)
        seg_rx = np.arange(seg_nx) - (seg_nx - 1)/2
        seg_ry = (seg_ny - 1)/2
        # loop through segments
        np.random.seed(seed)
        for i in range(seg_nx):
            seg_x =  seg_rx[i] * seg_d * np.cos(np.pi/6) 
            seg_y = -seg_ry[i] * seg_d
            for j in range(1, seg_ny[i]+1):
                # removes secondary and if any missing segment is present
                if (np.sqrt(seg_x**2 + seg_y**2) <= 4.01*seg_d) \
                        or ((seg_rx[i], j) in seg_missing):
                    seg_y += seg_d
                else:
                    # creates one hexagonal segment at x, y position in meters
                    segment = create_hexagon(nhr, seg_r, seg_y, seg_x, sampling)
                    # multiply by segment reflectivity and add to segments
                    seg_refl = np.random.normal(1, seg_rms)
                    segments += segment*seg_refl
                    seg_y += seg_d
        # need to transpose, due to the orientation of hexagons in create_hexagon
        segments = segments.T
        # resize to npupil, and add to pupil
        segments = resize_img(segments, npupil)
        pup *= segments

    return pup