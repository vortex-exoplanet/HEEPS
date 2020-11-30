import numpy as np

def polar_coord(npupil):
    
    xi,yi = cart_coord(npupil)
    ri = np.abs(xi + 1j*yi)
    ti = np.angle(xi + 1j*yi)

    return ri, ti

def cart_coord(npupil):

    dn = 2/npupil
    pup_range = np.arange(-1, 1, dn) + dn/2
    xi, yi = np.meshgrid(pup_range, pup_range)

    return xi, yi