import numpy as np
import scipy.optimize as opt
import scipy.special as spe
import matplotlib.pyplot as plt
from matplotlib import cm
from skimage.transform import resize
import warnings

def resize_img(img, new_size, preserve_range=True, mode='reflect',
        anti_aliasing=True):
    ''' Resize an image. Handles even and odd sizes.
    '''
    requirement = "new_size must be an int or a tuple/list of size 2."
    assert type(new_size) in [int, tuple, list], requirement
    if type(new_size) is int:
        new_size = (new_size, new_size)
    else:
        assert len(new_size) is 2, requirement
    assert  img.ndim in [2, 3], 'image must be a frame (2D) or a cube (3D)'
    if img.ndim == 3:
        new_size = (len(img), *new_size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # when anti_aliasing=False, and NANs
        img = np.float32(resize(np.float32(img), new_size, \
            preserve_range=preserve_range, mode=mode, anti_aliasing=anti_aliasing))
    
    return img

def pad_img(img, padded_size, pad_value=0):
    ''' Pad an img with a value (default is zero). Handles even and odd sizes.
    '''
    requirement = "padded_size must be an int or a tuple/list of size 2."
    assert type(padded_size) in [int, tuple, list], requirement
    if type(padded_size) is int:
        (x1, y1) = (padded_size, padded_size)
    else:
        assert len(padded_size) is 2, requirement
        (x1, y1) = padded_size
    (x2, y2) = img.shape
    # determine padding region
    assert not (x1<x2 or y1<y2), "padding region can't be smaller than image size."
    dx = int((x1 - x2)/2)
    dy = int((y1 - y2)/2)
    padx = (dx, dx) if (x1-x2)%2==0 else (dx+1, dx)
    pady = (dy, dy) if (y1-y2)%2==0 else (dy+1, dy)
    # pad image
    img = np.pad(img, [padx, pady], mode='constant', constant_values=pad_value)
    
    return img

def crop_img(img, new_size, margin=0):
    ''' Crop an img to a new size. Handles even and odd sizes. 
    Can add an optional margin of length 1, 2 (x,y) or 4 (x1,x2,y1,y2).
    '''
    requirement = "new_size must be an int or a tuple/list of size 2."
    assert type(new_size) in [int, tuple, list], requirement
    if type(new_size) is int:
        (x1, y1) = (new_size, new_size)
    else:
        assert len(new_size) is 2, requirement
        (x1, y1) = new_size
    (x2, y2) = img.shape
    # determine cropping region
    dx = int((x2 - x1)/2)
    dy = int((y2 - y1)/2)
    cropx = (dx, dx) if (x2-x1)%2==0 else (dx+1, dx)
    cropy = (dy, dy) if (y2-y1)%2==0 else (dy+1, dy)
    # check for margins
    requirement2 = "margin must be an int or a tuple/list of size 2 or 4."
    assert type(margin) in [int, tuple, list], requirement2
    if type(margin) is int:
        (mx1, mx2, my1, my2) = (margin, margin, margin, margin)
    elif len(margin) is 2:
        (mx1, mx2, my1, my2) = (margin[0], margin[0], margin[1], margin[1])
    else:
        assert len(margin) is 4, requirement2
        (mx1, mx2, my1, my2) = margin
    # crop image
    img = img[cropx[0]-mx1:-cropx[1]+mx2, cropy[0]-my1:-cropy[1]+my2]
    
    return img

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    ''' Model function. 2D Gaussian.
    '''
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
    
def oneD_Gaussian(x, amplitude, xo, sigma_x):
    ''' Model function. 1D Gaussian.
    '''
    xo = float(xo)
    #g = offset + amplitude*np.exp( ((x-xo)/(2*sigma_x))**2 )
    g = amplitude*np.exp( -((x-xo)/(np.sqrt(2)*sigma_x))**2 )
    #print(amplitude, xo, sigma_x)
    
    return g

def poly6(x, x0, a0, a1, a2, a3, a4, a5, a6):
    ''' Model function. Polynomial function up to 6th order.
    '''
    
    xx = x-x0
    y = a0 + a1*xx + a2*xx**2 + a3*xx**3 + a4*xx**4 + a5*xx**5 + a6*xx**6
    
    return y

def poly6odd(x, x0, a0, a2, a4, a6):
    ''' Model function. Polynomial function up to 6th order.
    '''
    xx = x-x0
    y = a0 + a2*xx**2 + a4*xx**4 + a6*xx**6
    return y

def twoD_Airy(xy, amplitude, xo, yo, F):
    ''' Model function. 2D Airy.
    '''
    
    (x, y) = xy
    r = np.sqrt((x-xo)**2+(y-yo)**2)*F
    
    nx=r.shape[1]
    ny=r.shape[0]
    maxmap=np.where(r==0, np.ones((ny,nx)), np.zeros((ny,nx)))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.unravel_index(maxmap.argmax(), maxmap.shape)
        r[indmax]=1.
    elif nbmax > 1:
        print('ERROR in twoD_Airy: several nulls')
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1:
        Airy[indmax]=amplitude
    
    return Airy.ravel()

def oneD_Airy(x, amplitude, xo, F):
    ''' Model function. 1D Airy.
    '''
    
    r=(x-xo)*F
    nx=x.shape[0]
    
    maxmap=np.where(x==0, np.ones(nx), np.zeros(nx))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.argmax(maxmap)
        r[indmax]=1.
    elif nbmax > 1:
        print('ERROR in oneD_Airy: several nulls')
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1:
        Airy[indmax]=amplitude
    
    return Airy

def oneD_Airy_log(x, amplitude, xo, F):
    ''' Model function. 1D log10(Airy).
    '''    
    
    r=(x-xo)*F
    nx=x.shape[0]
    
    maxmap=np.where(r==0, np.ones(nx), np.zeros(nx))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.argmax(maxmap)
        r[indmax]=1.
    elif nbmax > 1:
        print('ERROR in oneD_Airy: several nulls')
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1 :   
        Airy[indmax]=amplitude
    
    return np.log10(Airy)

def fit_gauss_2D(img):
    ''' Fits a 2D Gaussian pattern on the image.
        
        Returns the best fit parameters of the Gaussian shape.
        
        See twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''
    
    nx=img.shape[1]
    ny=img.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    x, y = np.meshgrid(x, y)
    xy = (x, y)
    
    init_xmax=np.unravel_index(img.argmax(), img.shape)[1]
    init_ymax=np.unravel_index(img.argmax(), img.shape)[0]
    initial_guess = (img.max(), init_xmax, init_ymax, 5, 5, 0, 0)
    popt, pcov = opt.curve_fit(twoD_Gaussian, xy, img.ravel(), 
                               p0=initial_guess)
    
    return popt

def fit_gauss_1D(y, x):
    ''' Fits a 1D Gaussian curve.
        
        Returns the best fit parameters of the Gaussian shape.
        
        See oneD_Gaussian(x, amplitude, xo, sigma_x, offset)
    '''
    
    #nx=y.shape[0]
    #x = np.linspace(0, nx-1, nx)
    init_xmax=x[y.argmax()]
    
    initial_guess = (y.max(), init_xmax, (x[-1]-x[0])/4.)
    popt, pcov = opt.curve_fit(oneD_Gaussian, x, y, p0=initial_guess)
    
    return popt

def fit_airy_2D(img, disp=0):
    ''' Fits a 2D Airy pattern on the image.
        
        Returns the best fit parameters of the Airy pattern.
        
        See twoD_Airy(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''
    
    nx=img.shape[1]
    ny=img.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    x, y = np.meshgrid(x, y)
    xy = (x, y)
    
    init_xmax=np.unravel_index(img.argmax(), img.shape)[1]
    init_ymax=np.unravel_index(img.argmax(), img.shape)[0]
    initial_guess = (img.max(), init_xmax, init_ymax, .4)
    #initial_guess = (img[init_xmax, init_ymax]  , init_xmax, init_ymax, .6)
    #plt.figure(27)
    #plt.imshow(twoD_Airy(xy, img[init_xmax, init_ymax]  , init_xmax, \
    #        init_ymax, .6).reshape(nx,ny))
    
    popt, pcov = opt.curve_fit(twoD_Airy, xy, img.ravel(), p0=initial_guess)
    
    if disp != 0:
        data_fitted = twoD_Airy(xy, *popt)
        #offset=np.min(img)
        plt.figure(disp)
        plt.clf()
        plt.subplot(1,3,1)
        plt.imshow(data_fitted.reshape(nx,ny), interpolation='none', \
                cmap=cm.Greys_r)
        plt.colorbar()
        plt.subplot(1,3,2)
        plt.plot(img[popt[1],:])
        plt.plot(data_fitted.reshape(nx,ny)[popt[1],:], 'r--')
        plt.yscale('log')
        plt.subplot(1,3,3)
        plt.plot(img[:,popt[2]])
        plt.plot(data_fitted.reshape(nx,ny)[:,popt[2]], 'r--')
        plt.yscale('log')
        
        plt.figure(disp+1)
        plt.clf()
        plt.subplot(121)
        plt.imshow(data_fitted.reshape(nx,ny), interpolation='none', \
                cmap=cm.Greys_r)
        plt.colorbar()
        P_fit=get_radial_profile(data_fitted.reshape(nx,ny), (popt[1], \
                popt[2]), 1, disp=10)
        P_mes=get_radial_profile(img, (popt[1], popt[2]), 1, disp=0)
        plt.subplot(122)
        plt.plot(P_mes, 'b')
        plt.plot(P_fit, 'r--')
        plt.yscale('log')
        
        print('\n--- Airy disk fit results ---')
        print('Amplitude ='+str(popt[0]))
        print('Position of the maximum: \nxo='+str(popt[1])+' nyo='+str(popt[2]))
        print('F factor='+str(popt[3]))
        print('-----------------------------')
    
    return popt

def fit_airy_1Dlog(Y, disp=0, initial_guess=[0.,0.,0.]):
    "fit one D   "
    
    nx=Y.shape[0]
    x=np.linspace(0,nx-1,nx)
    
#    minval=np.min(Y)
#    Y-=minval
#    ampl=np.max(Y)
#    Y=Y/ampl
    
    Ylog=np.log10(Y)
    
    if np.sum(initial_guess) == 0. :
        initial_guess=(np.max(Y), np.argmax(Y), .5)
    
    popt, pcov = opt.curve_fit(oneD_Airy_log, x, Ylog, p0=initial_guess, \
            sigma=1./Y**2)
    
    if disp != 0:
        data_fitted=oneD_Airy_log(x, *popt)
        plt.figure(disp)
        plt.clf()
        plt.plot(Ylog, 'b')
        plt.plot(x, data_fitted, 'r--')
    
    return popt

def fit_airy_1D(Y, disp=0, initial_guess=[0.,0.,0.]):
    "fit one D   "
    
    nx=Y.shape[0]
    x=np.linspace(0,nx-1,nx)
    
#    minval=np.min(Y)
#    Y-=minval
#    ampl=np.max(Y)
#    Y=Y/ampl
    
    if np.sum(initial_guess) == 0. :
        initial_guess=(np.max(Y), np.argmax(Y), .5)
    
    popt, pcov = opt.curve_fit(oneD_Airy, x, Y, p0=initial_guess, sigma=1./Y)
    
    if disp != 0:
        data_fitted=oneD_Airy(x, *popt)
        plt.figure(disp)
        plt.clf()
        plt.plot(Y, 'b')
        plt.yscale('log')
        plt.plot(x, data_fitted, 'r--')
    
    return popt
    
def get_r_dist(nx,ny,xo,yo):
    ''' Returns the array of dimensions (nx,ny) with values corresponding the
        distance from the center (xo,yo). 
    '''
    
    x = np.linspace(0, nx, nx)-xo-1
    y = np.linspace(0, ny, ny)-yo-1
    x, y = np.meshgrid(x, y)
    
    return np.sqrt(x**2+y**2)

def get_radial_profile(img, xoyo, nbin, disp=0):
    ''' Computes the mean radial profile of the image.
    
        img:
            2D image.
        xoyo:
            (xo,yo) center for the annuli.
        nbin:
            width of the annuli in pixels
        disp:
            optional key word for displaying the images.
            Its value will serve as the window number that will be created.
    '''
    
    (xo, yo) = xoyo
    (nx, ny) = img.shape
    r = get_r_dist(nx, ny, xo, yo)
    
    r_max = np.max(r) # radius of the image
    r_max = np.max(r[xo,:])
    
    npts=int(r_max/nbin)
    O=np.ones((nx,ny))
    Z=np.zeros((nx,ny))
    Profile=np.zeros(npts)
    
    if disp != 0:
        plt.figure(disp)
        plt.clf()
        plt.subplot(121)
        plt.plot(xo,yo, 'xw')
        plt.title('PSF')
        plt.title('Averaged radial profile')
        plt.xlabel('Distance from center (pixels)')
        plt.imshow(img, interpolation='none')
        plt.colorbar()
        
        val_min=np.min(img)
        val_max=np.max(img)
        
        for k in range(0,npts-1,1):
            M=np.where(r>nbin*k, O, Z)*np.where(r<nbin*(k+1), O, Z)
            Profile[k]=np.sum(img*M)/np.sum(M)
            
            plt.figure(disp)
            plt.subplot(121)
            plt.imshow(img, interpolation='none')
            plt.pause(.005)
            plt.imshow(img*M, interpolation='none', vmin=val_min, vmax=val_max)
            plt.pause(.005)
            plt.subplot(122)
            plt.plot(Profile, 'rx')
            plt.yscale('log')
        
        plt.plot(Profile, 'r')
        plt.subplot(121)
        plt.imshow(img, interpolation='none', vmin=val_min, vmax=val_max)
    
    for k in range(0,npts-1,1):
        M=np.where(r>nbin*k, O, Z)*np.where(r<nbin*(k+1), O, Z)
        Profile[k]=np.sum(img*M)/np.sum(M)   
    
    return Profile

def adjust_bckgr_level(img, xo, yo, R=0, disp=0):
    ''' Computes the median/mean background level of the image outside a 
        given radius.
        
        img:
            2D image
        (xo,yo):
            center of the PSF
        R: 
            radius of the circular zone to exclude.
        disp:
            optional keyword for displaying images.
    '''
    
    (nx,ny)=img.shape
    r=get_r_dist(nx,ny,xo,yo)
    
    M=np.double(img[np.where(r>R)])   
    
    bckgr_med=np.median(M)
    bckgr_mean=np.mean(M)
    
    
    print('\n----- Checking background level ------')
    print('Background level = '+str(bckgr_med)+' [median]')
    print('                   '+str(bckgr_mean)+' [mean]')
    
    if disp != 0:
        plt.figure(disp)
        plt.clf()
        Pattern=np.where(r>R,np.ones((nx,ny)),np.zeros((nx,ny)))
        plt.imshow(Pattern*img, interpolation='none')
        plt.colorbar()
        plt.title('Background median on external area = '+str(bckgr_med))
    
    if bckgr_med != 0 :
        img=img-bckgr_med
        print(' -> background adjusted to 0 median')
    
    print('------------------------------------- ')
        
    return (img, bckgr_med, bckgr_mean)
