import numpy as np

def create_hexagon(ngrid, hexrad, xc, yc, sampling, mag=11):
    """Returns an image containing a hexagon with antialiased edges.
            at postion xc and yc in meters
       
    Input
    ----------    
    ngrid : int
        grid size
        
    hexrad : float
        Distance from any vertex to the center of the hexagon
    
    xc, yc : float
        center of the hexagon in meters
   
    sampling: flaot
        meters/pixel
    
    mag: int
        sampling factor to create anti-aliased image
    
    Returns
    -------
    image : numpy ndarray
        Returns an image containing a hexagon.
    """
    # take a larger grid, and crop it at the end
    ngrid *= 2
    
    image = np.zeros([ngrid, ngrid], dtype = np.float64)
    xcpix = xc / sampling + ngrid/2 - 1
    ycpix = yc / sampling + ngrid/2 - 1
    radpix = hexrad/sampling
    
    # angles for verices of hexagon
    angles = np.array([90., 30., -30., -90., 210, 150])
    
    xv = np.cos(np.deg2rad(angles))
    yv = np.sin(np.deg2rad(angles))
    
    xp = xv * radpix + xcpix 
    yp = yv * radpix + ycpix
        
    yindex = np.round([yp[3], yp[4], yp[5], yp[0]+1]).astype(int)
    
    slopes = 1 / np.tan(np.deg2rad([30, 90, 150]))
    
    for l in range(len(yindex)-1):
        for ypix in range(yindex[l], yindex[l+1]):
            mag2 = mag
            for ysub in range(0, mag):
                y = ypix - 0.5 + (0.5 + ysub)/mag # this breaks a pixel into mag steps
                
                if y < yp[3]:
                    continue # if true skips rest of the code, to exclude values y<yp[left]
                if y > np.max(yp[0]):
                    break # this is to break the magnification loop when y reaches max(yp)
                   
                if y<yp[4]:
                    ind = 0
                if y < yp[5] and y >= yp[4]:
                    ind = 1
                if y < yp[0] and y >= yp[5]:
                    ind = 2
                xleft = -slopes[ind] * (y-yp[3+ind]) + xp[3+ind]
                xright = slopes[ind] * (y-yp[3-ind]) + xp[3-ind]
                xleftpix = int(np.round(xleft))
                xrightpix = int(np.round(xright))
                
                if xleftpix != xrightpix:
                    image[ypix,xleftpix] = image[ypix,xleftpix] + mag * ((xleftpix + 0.5) - xleft)
                    image[ypix,xrightpix] = image[ypix,xrightpix] + mag * (xright - (xrightpix - 0.5))
                    
                    if (xrightpix - xleftpix) > 1:
                        imin = xleftpix+1
                        imax = xrightpix
                        image[ypix,imin:imax] +=  mag
                else:
                    image[ypix,xleftpix] = image[ypix,xleftpix] + mag * (xright - xleft)
    
    image = image / float(mag)**2
    
    # crop final image
    crop = int(ngrid/4)
    image = image[crop:3*crop, crop:3*crop]
        
    return image

