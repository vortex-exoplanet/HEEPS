import numpy as np

def ri_ti(npupil):
    
    dn = 2/npupil
    pup_range = np.arange(-1, 1, dn) + dn/2
    xi,yi = np.meshgrid(pup_range, pup_range)
    ri = np.abs(xi + 1j*yi)
    ti = np.angle(xi + 1j*yi)

    return ri, ti