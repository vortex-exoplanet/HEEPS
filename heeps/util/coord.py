import numpy as np

def polar_coord(npupil):
    
    x, y = cart_coord(npupil)
    r = np.abs(x + 1j*y)
    t = np.angle(x + 1j*y)

    return r, t

def cart_coord(npupil):

    dn = 2/npupil
    pup_range = np.arange(-1, 1, dn) + dn/2
    x, y = np.meshgrid(pup_range, pup_range)

    return x, y

def disk_coord(radius, nr=4):
    
    dr = radius/(nr - 0.5)
    x, y = [0], [0] # center
    for r in np.arange(radius/dr)*dr:
        npts = int(round(2*np.pi*r/dr))
        for t in np.arange(npts)*2*np.pi/npts:
            z = r*np.exp(1j*t)
            x.append(np.real(z))
            y.append(np.imag(z))
    
    return x, y