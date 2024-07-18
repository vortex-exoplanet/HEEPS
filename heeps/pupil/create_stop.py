from heeps.pupil import create_pupil
from heeps.util.coord import cart_coord, polar_coord
from heeps.util.img_processing import resize_img
import numpy as np
import proper

def create_stop(d_ext, d_int, dRext, dRint, dRspi, nhr=1023, npupil=285, ngrid=1024,
        pupil_img_size=40, diam_nominal=38, spi_width=0.54, seg_width=1.45,
        spi_angles=[0,60,120,180,240,300], AP_width=0, AP_angles=[0], 
        AP_center=0.5, misalign_x=0, misalign_y=0, circ_ext=True, 
        circ_int=True, **conf):
    '''
    Margins are calculated wrt nominal diameter, and applied (added/subtracted) 
    to the external/internal diameters, and spider width.
    '''
    # save/store PROPER ngrid value
    ntmp = proper.n
    # nhr must be >= npupil
    if nhr < npupil:
        nhr = 10*npupil-1 # odd
    # misalignments
    dx = misalign_x*diam_nominal/pupil_img_size
    dy = misalign_y*diam_nominal/pupil_img_size
    # create spider stop
    conf_spi = dict(
        npupil = nhr,               # high res
        pupil_img_size = pupil_img_size, 
        diam_ext = 2*pupil_img_size,# no circular aperture
        diam_int = 0,               # no central obscuration
        seg_width = 0,              # no segments
        spi_angles = spi_angles,
        spi_width = spi_width + dRspi*diam_nominal,
        dx = dx,
        dy = dy)
    mask_spi = create_pupil(**conf_spi)
    # create thicker spiders for asym pupil (AP)
    if AP_width > 0:
        conf_spi.update(
            spi_angles = AP_angles,
            spi_width = AP_width,
            spi_norm_center = AP_center)
        mask_spi *= create_pupil(**conf_spi)
    # calculate dodecagonal (outer) diameter
    if circ_ext == False:
        alpha = np.arcsin(seg_width/d_ext)
        d_ext *= np.cos(alpha)
    # calculate dodecagonal (inner) diameter
    if circ_int == False:
        beta = np.pi/6 - np.arcsin(seg_width*np.sin(np.pi/6)/d_int)
        d_int *= np.cos(beta)
    # create outer and inner stops
    r_ext = (d_ext - dRext*diam_nominal) / pupil_img_size
    r_int = (d_int + dRint*diam_nominal) / pupil_img_size
    # create mask
    r, t = polar_coord(nhr, dx=dx, dy=dy)
    mask_ext = (r < r_ext) if circ_ext == True \
                         else 1 - dodecagon(r_ext, nhr, dx=dx, dy=dy)
    mask_int = (r > r_int) if circ_int == True \
                         else 1 - hexagon(r_int, nhr, dx=dx, dy=dy)
    # resize
    mask = resize_img(mask_ext*mask_int*mask_spi, npupil)
    # reload PROPER saved ngrid value
    proper.n = ntmp
    return mask

def dodecagon(r_ext, npupil, dx=0, dy=0):
    x, y = cart_coord(npupil, dx=dx, dy=dy)
    M1 = mask_angle(x, y, r_ext, 0)
    M2 = mask_angle(x, y, r_ext, 30)
    M3 = mask_angle(x, y, r_ext, 60)
    M4 = mask_angle(x, y, r_ext, 90)
    return np.uint8(M1 + M2 + M3 + M4)

def hexagon(r_int, npupil, dx=0, dy=0):
    x, y = cart_coord(npupil, dx=dx, dy=dy)
    M5 = mask_angle(x, y, r_int, 0)
    M6 = mask_angle(x, y, r_int, 60)
    return 1 - np.uint8(M5 + M6)

def mask_angle(x, y, r, theta):
    P = (r*np.cos(np.deg2rad(theta)), r*np.sin(np.deg2rad(theta)))
    a = np.tan(np.deg2rad(theta - 90))
    b = P[1] - a*P[0]
    mask = (y > a*x + b) + (y > -a*x + b) + (y < a*x - b) + (y < -a*x - b)
    return mask