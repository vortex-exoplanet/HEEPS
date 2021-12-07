import numpy as np
import astropy.units as u

def cart_coord(npupil, dx=0, dy=0):
    ''' dx, dy: misalignment percentage 
    '''
    dn = 2/npupil
    pup_range = np.arange(-1, 1, dn) + dn/2
    x, y = np.meshgrid(pup_range - 2*dx, pup_range - 2*dy)
    return np.float32(x), np.float32(y)

def polar_coord(npupil, dx=0, dy=0):
    x, y = cart_coord(npupil, dx=dx, dy=dy)
    r = np.abs(x + 1j*y)
    t = np.angle(x + 1j*y)
    return np.float32(r), np.float32(t)

def disk_coord(radius, nr=4):
    dr = radius/(nr - 0.5)
    x, y = [0], [0] # center
    for r in np.arange(radius/dr)*dr:
        npts = int(round(2*np.pi*r/dr))
        for t in np.arange(npts)*2*np.pi/npts:
            z = r*np.exp(1j*t)
            x.append(np.real(z))
            y.append(np.imag(z))
    return np.array([x, y]).T

def mas2rms(xy, diam):
    '''
    translate the tip/tilt to RMS phase errors
    tilt in mas: RMS = x/2 = (tilt*(D/2))/2 = tilt*diam/4
    tilt in lam/D -> tilt*diam/4 *(lam/diam) = tilt*lam/4
    '''
    return np.array(xy)*u.mas.to('rad')*diam/4

def rms2mas(xy, diam):
    return np.array(xy)*u.rad.to('mas')*4/diam