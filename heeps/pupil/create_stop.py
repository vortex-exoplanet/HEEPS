from .create_hexagon import create_hexagon
from heeps.util.coord import cart_coord
from heeps.util.img_processing import resize_img
import proper
import numpy as np


def create_stop(ngrid=2**10, npupil=285, pupil_img_size=40, diam_nominal=38,
        diam_ext=37, diam_int=11, spi_angles=[0,60,120,180,240,300],
        spi_width=0.54, AP_angles=[], AP_width=0, AP_center=1, AP_length=2,
        dRext=0, dRint=0, dRspi=0, force_sym=False, dx=0, dy=0,
        circ_ext=True, circ_int=True, add_seg=False,
        seg_ny=[10,13,16,19,22,23,24,25,26,27,28,29,30,31,30,31,
        30,31,30,31,30,31,30,29,28,27,26,25,24,23,22,19,16,13,10],
        seg_missing=[], seg_rint=4.01, seg_width=1.45, seg_gap=0.004, seg_ptv=0,
        seed=123456, verbose=False, **conf):
    
    ''' Create a pupil.
    
    Args:
        ngrid: int
            number of pixels of the wavefront array
        npupil: int
            number of pixels of the pupil
        pupil_img_size: float
            pupil image (for PROPER) in m
        diam_nominal: float
            nominal diameter in m
        diam_ext: float
            outer circular aperture in m
        diam_int: float
            central obscuration in m
        spi_angles: list of float
            refular spider angles in deg
        spi_width: float
            regular spider width in m
        AP_angles: list of float
            asymetric pupil spider angles in deg
        AP_width: float
            asymetric pupil spider width in m (AP_width > spi_width)
        AP_center: float
            asymetric pupil spider normalized center wrt nominal radius
        AP_length: float
            asymetric pupil spider normalized length wrt nominal radius
        dRext, dRint, dRspi: floats
            margins for diam_ext, diam_int, spi_width, wrt diam_nominal
        force_sym : bool
            if True and dRspi is not 0, force the spiders to have all the same dimension,
            max(spi_width, AP_width). Specially useful to force symmetric Lyot stop.
            Default: False
        dx, dy: floats
            percentage of misalignment in x and y
        circ_ext, circ_int: bools
            True if circular aperture/obscuration
        add_seg: bool
            True if adding segments
        seg_ny: list of int
            number of hexagonal segments per column (from left to right)
        seg_missing: list of tupples
            coordinates of missing segments
        seg_rint: float
            number of black segments accross the central obstruction radius
        seg_width: float
            segment width in m
        seg_gap: float
            gap between segments in m
        seg_ptv: float
            ptv uniform segment intensity loss
        seed: int
            seed used to initialize the random number generator

    '''

    if verbose is True:
        print("dRext=%s, dRint=%s, dRspi=%s, circ_ext=%s, circ_int=%s, AP_angles=%s, add_seg=%s"\
            %(dRext, dRint, dRspi, circ_ext, circ_int, AP_angles, add_seg))

    # store PROPER ngrid value
    ntmp = proper.n

    # create a high res pupil with PROPER
    ngrid += ngrid % 2 # ngrid must be even
    while ngrid < npupil:
        ngrid *= 2
    # one row and column will be removed later: pupil_img_size => (ngrid-1)
    wf_tmp = proper.prop_begin(1, 1, ngrid, (diam_nominal/pupil_img_size) * (ngrid-1)/ngrid) 

    # external aperture
    if diam_ext > 0:
        if circ_ext == True:
            proper.prop_circular_aperture(wf_tmp, (diam_ext/2) / (diam_nominal/2) - dRext, 
                                          dx, dy, NORM=True) # NORM: nominal radius = 1
        else: # dodecagonal aperture for METIS
            alpha = np.arcsin(seg_width/diam_ext)
            diam_ext *= np.cos(alpha) # small corrrection for diam_ext
            r_ext = (diam_ext - dRext*diam_nominal) / pupil_img_size
            mask_ext =  1 - dodecagon(r_ext, ngrid-1, dx=dx, dy=dy)
            # add one extra row and column of zeros
            mask_ext = np.hstack([np.zeros((ngrid, 1)), np.vstack([np.zeros((1, ngrid-1)), mask_ext])])
            proper.prop_multiply(wf_tmp, mask_ext)

    # central obscuration
    if diam_int > 0:
        if circ_int == True:
            proper.prop_circular_obscuration(wf_tmp, (diam_int/2) / (diam_nominal/2) + dRint, 
                                             dx, dy, NORM=True) # NORM: nominal radius = 1
        else: # hexagonal obscuration for METIS
            beta = np.pi/6 - np.arcsin(seg_width*np.sin(np.pi/6)/diam_int)
            diam_int *= np.cos(beta) # small corrrection for diam_int
            r_int = (diam_int + dRint*diam_nominal) / pupil_img_size
            mask_int =  1 - hexagon(r_int, ngrid-1, dx=dx, dy=dy)
            # add one extra row and column of zeros
            mask_int = np.hstack([np.zeros((ngrid, 1)), np.vstack([np.zeros((1, ngrid-1)), mask_int])])
            proper.prop_multiply(wf_tmp, mask_int)

    if (dRspi != 0) and force_sym:
        # Lyot stop oversize handling in case we want to force the same spiders sizes.
        print('[Warning] Forcing spiders to have all the same sizes.')
        new_spi_width = np.max([spi_width, AP_width])
        new_spi_angles = set(spi_angles)     # Convert the original list to a set to remove duplicates
        new_spi_angles.update(AP_angles)    # Add elements from the new list to the set
        new_spi_angles = list(new_spi_angles)    # Convert the set back to a list
        if len(new_spi_angles) > 0:
            for angle_deg in new_spi_angles:
                angle_rad = np.radians(angle_deg)
                spi_length = pupil_img_size/diam_nominal
                proper.prop_rectangular_obscuration(wf_tmp, 2 * (new_spi_width/diam_nominal + dRspi), spi_length,
                    np.sin(angle_rad)*spi_length/2 + dx, -np.cos(angle_rad)*spi_length/2 + dy,
                    ROTATION=angle_deg, NORM=True) # NORM: nominal radius = 1
    else:
        # regular spiders
        if len(spi_angles) > 0:
            for angle_deg in spi_angles:
                angle_rad = np.radians(angle_deg)
                spi_length = pupil_img_size/diam_nominal
                proper.prop_rectangular_obscuration(wf_tmp, 2 * (spi_width/diam_nominal + dRspi), spi_length,
                    np.sin(angle_rad)*spi_length/2 + dx, -np.cos(angle_rad)*spi_length/2 + dy,
                    ROTATION=angle_deg, NORM=True) # NORM: nominal radius = 1
        
        # asymetric pupil spiders
        if len(AP_angles) > 0:
            for angle_deg in AP_angles:
                angle_rad = np.radians(angle_deg)
                proper.prop_rectangular_obscuration(wf_tmp, 2 * (AP_width/diam_nominal + dRspi), AP_length,
                    np.sin(angle_rad)*AP_center + dx, -np.cos(angle_rad)*AP_center + dy,
                    ROTATION=angle_deg, NORM=True) # NORM: nominal radius = 1

    # crop the pupil to pupil_img_size => (ngrid-1)
    pup = proper.prop_get_amplitude(wf_tmp)[1:,1:]
    # resize to npupil
    pup = resize_img(pup, npupil)

    # add segments
    if add_seg is True:
        segments = np.zeros((ngrid, ngrid))
        # sampling in meters/pixel
        sampling = pupil_img_size/ngrid
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
                if (np.sqrt(seg_x**2 + seg_y**2) <= seg_rint*seg_d) \
                        or ((seg_rx[i], j) in seg_missing):
                    pass
                else:
                    # creates one hexagonal segment at x, y position in meters
                    segment = create_hexagon(ngrid, seg_r, seg_y, seg_x, sampling)
                    # calculate segment reflectivities in amplitude (seg_ptv=intensity)
                    seg_refl = np.random.uniform(1-seg_ptv, 1)**0.5
                    # multiply, then add segment to segments
                    segments += segment*seg_refl
                # go to next segment of the column
                seg_y += seg_d
        # need to transpose, due to the orientation of hexagons in create_hexagon
        segments = segments.T
        # resize to npupil, and add to pupil
        segments = resize_img(segments, npupil)
        pup *= segments

    # reload PROPER saved ngrid value
    proper.n = ntmp

    return pup

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