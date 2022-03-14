from .create_pupil import create_pupil
from heeps.util.coord import cart_coord, polar_coord
import numpy as np


def create_ls(d_ext, d_int, npupil=285, pupil_img_size=40, diam_nominal=38, 
        spi_width=0.5, ls_dRext=0, ls_dRint=0, ls_dRspi=0, ls_misalign_x=0, 
        ls_misalign_y=0, circ_ext=True, circ_int=True, **conf):
    # misalignments
    dx = ls_misalign_x*diam_nominal/pupil_img_size
    dy = ls_misalign_y*diam_nominal/pupil_img_size
    # create spider stop
    conf = dict(
        npupil=npupil,
        pupil_img_size=pupil_img_size, 
        diam_ext=2*pupil_img_size,  # no circular aperture
        diam_int=0,                 # no central obscuration
        seg_width=0,                # no segments
        spi_width=spi_width + ls_dRspi*diam_nominal,
        dx=dx,
        dy=dy)
    mask_spi = create_pupil(**conf)
    # create outer and inner stops
    r_ext = (d_ext - ls_dRext*diam_nominal) / pupil_img_size
    r_int = (d_int + ls_dRint*diam_nominal) / pupil_img_size
    # create mask
    r, t = polar_coord(npupil, dx=dx, dy=dy)
    mask_ext = (r < r_ext) if circ_ext == True \
                         else ~dodecagon(r_ext, npupil, dx=dx, dy=dy)
    mask_int = (r > r_int) if circ_int == True \
                         else ~hexagon(r_int, npupil, dx=dx, dy=dy)
    # lyot stop
    return mask_ext * mask_int * mask_spi

def dodecagon(r_ext, npupil, dx=0, dy=0):
    x, y = cart_coord(npupil, dx=dx, dy=dy)
    M1 = mask_angle(x, y, r_ext, 0 , 90)
    M2 = mask_angle(x, y, r_ext, 30, 90)
    M3 = mask_angle(x, y, r_ext, 60, 90)
    M4 = mask_angle(x, y, r_ext, 90, 90)
    return M1 + M2 + M3 + M4

def hexagon(r_int, npupil, dx=0, dy=0):
    x, y = cart_coord(npupil, dx=dx, dy=dy)
    M5 = mask_angle(x, y, r_int, 30, 60)
    M6 = mask_angle(x, y, r_int, 30, 120)
    return ~(M5 + M6)

def mask_angle(x, y, r, theta, rot):
    P = (r*np.cos(np.deg2rad(theta)), r*np.sin(np.deg2rad(theta)))
    a = np.tan(np.deg2rad(theta - rot))
    b = P[1] - a*P[0]
    mask = (y > a*x + b) + (y > -a*x + b) + (y < a*x - b) + (y < -a*x - b)
    return mask