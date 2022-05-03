from heeps.util.img_processing import resize_cube, get_rms
from heeps.util.multiCPU import multiCPU
import astropy.convolution as astroconv
from astropy.io import fits
import numpy as np
import warnings
import proper
from copy import deepcopy
from scipy.signal import savgol_filter

def conv_kernel(Npup, cpp, HR=2**11):
    ker_range = np.arange(-HR, HR)/HR
    Xs,Ys = np.meshgrid(ker_range, ker_range)       # high res kernel XY grid
    Rs = np.abs(Xs + 1j*Ys)                         # high res kernel radii
    kernel = np.ones((2*HR, 2*HR))
    kernel[Rs > 1] = 0
    # kernel must have odd dimensions
    nkernel = int(Npup/cpp)
    nkernel = nkernel+1 if nkernel%2 == 0 else nkernel
    # resize kernel
    kernel = resize_cube(kernel, nkernel)
    kernel /= np.sum(kernel)                        # need to normalize the kernel
    return kernel

def spatial(allSF, kernel, npupil=None, norm=False, verbose=False):
    # mask with nans
    mask_nan = np.isnan(allSF) + (allSF == 0)
    allSF[mask_nan] = np.nan
    # get low and high spatial frequencies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # NANs
        LSF = astroconv.convolve(allSF, kernel, boundary='extend')
        LSF[mask_nan] = np.nan
        HSF = allSF - LSF
    # print rms
    if verbose is True:
        print('rms(all SF) = %.2f nm'%(np.nanstd(allSF)))
        print('rms(LSF) = %.2f nm'%(np.nanstd(LSF)))
        print('rms(HSF) = %.2f nm'%(np.nanstd(HSF)))
    # normalize
    if norm is True:
        allSF /= np.nanstd(allSF)
        LSF /= np.nanstd(LSF)
        HSF /= np.nanstd(HSF)
        allSF -= np.nanmean(allSF)
        LSF -= np.nanmean(LSF)
        HSF -= np.nanmean(HSF)
    # remove nans
    allSF = np.nan_to_num(allSF)
    LSF = np.nan_to_num(LSF)
    HSF = np.nan_to_num(HSF)
    # resize outputs
    if npupil is not None:
        allSF = resize_cube(allSF, npupil)
        LSF = resize_cube(LSF, npupil)
        HSF = resize_cube(HSF, npupil)
    return allSF, LSF, HSF

def temporal(t_max, dt, fc1, fc2, seed=123456):
    '''Parseval's Theorem: sum of squares in the time domain equals sum 
    of squares in the frequency domain (or total energy in the time domain 
    equals total energy in the frequency domain)'''
    np.random.seed(seed)
    # Time domain (time series)
    N = int(t_max/dt)                   # number of samples in the time series (e.g. 12000)
    ts = np.arange(N + 1) * dt          # time steps
    # Frequency domain (power spectral density)
    df = 1 / t_max                      # sampling in Hz
    M = int(N/2)                        # number of samples in the PSD
    fs = np.arange(M + 1) * df          # frequency steps
    f_max = M * df                      # maximum frequency (e.g. 1.66 Hz)
    assert 0 <= fc1 <= f_max, "lower cutoff is out of range."
    assert 0 <= fc2 <= f_max, "upper cutoff is out of range."
    rms = 1                            # rms value
    # Calculate the Inverse Fourier transform (time series)
    G = N * np.sqrt(df/2)                       # norm factor
    Lrect = fc2 - fc1                           # length of the rectangular function
    Mrect = (fs >= fc1)*(fs <= fc2)             # mask of the rectangular function
    Phis = 2*np.pi*np.random.random(len(fs))    # spectral phase
    Amps = rms*np.sqrt(1/Lrect)*Mrect           # spectral amplitude
    def complex_spectrum(Phis, Amps):
        # cf. https://stackoverflow.com/questions/9062387/ifft-of-symmetric-spectrum
        Yf =  Amps * np.exp(1j*Phis)
        Yf_left = np.append(0, Yf[:-1])
        Yf_right = np.flip(np.conj(Yf[:-1]))
        return np.append(Yf_left, Yf_right)
    Yf = complex_spectrum(Phis, Amps)           # complex spectrum
    ys = np.fft.ifft(Yf) * G                    # time series
    # select N real scaling factors, and apply rms correction
    ys1 = np.real(ys[:N])
    rms_corr = rms**2/np.std(ys1)
    Amps = rms_corr*np.sqrt(1/Lrect)*Mrect
    Yf = complex_spectrum(Phis, Amps)
    ys = np.fft.ifft(Yf) * G
    ys1 = np.real(ys[:N])
    return ys1

def fit_zer(pup, rad, nzer, cube):
    # pup should NOT have NANs for prop_fit_zernikes
    return proper.prop_fit_zernikes(cube, pup, rad, nzer, eps=0, fit=True)

def get_zernike(cube_name, pup, nzer):
    zpols_name = cube_name[:-5] + '_zpols_%s.fits'%nzer
    try:
        zpols = fits.getdata(zpols_name)
        print('getdata ' + zpols_name)
    except FileNotFoundError:
        print('writeto ' + zpols_name)
        cube = fits.getdata(cube_name)
        nimg = cube.shape[-1]
        pup = resize_cube(pup, nimg)
        zpols = multiCPU(fit_zer, posargs=[pup, nimg/2, nzer], 
                posvars=[cube], case='get zpols')
        fits.writeto(zpols_name, np.float32(zpols))
    return zpols

def remove_zernike(wf, pup, allSF, zpols):
    proper.prop_add_phase(wf, allSF)
    LSF = proper.prop_zernikes(wf, np.arange(len(zpols)) + 1, zpols).astype('float32')
    LSF[pup==0] = 0
    HSF = allSF.astype('float32') - LSF    
    return LSF, HSF

def psd_spatial_zernike(cube_name, pup, zpols, nzer, ncube):
    spsd_name = cube_name[:-5] + '_%s' + '_%s.fits'%nzer
    try:
        spsd = fits.getdata(spsd_name%'spsd')
        print('getdata ' + spsd_name%'spsd')
    except FileNotFoundError:
        print('writeto ' + spsd_name%'spsd')
        cube = fits.getdata(cube_name)[:ncube]
        nimg = cube.shape[-1]
        pup = resize_cube(pup, nimg)
        zpols = zpols[:ncube,:nzer]
        wf = proper.prop_begin(1, 1, nimg, 1) # initial wavefront
        LSFs = np.empty((nzer, ncube, nimg, nimg))
        HSFs = np.empty((nzer, ncube, nimg, nimg))
        HSFs_rms = []
        for z in np.arange(nzer) + 1:
            verbose = True if z == 1 else False
            LSF, HSF = multiCPU(remove_zernike, posargs=[deepcopy(wf), pup],
                        posvars=[cube, zpols[:,:z]], case='remove zernike', 
                        multi_out=True, verbose=verbose)
            LSFs[z-1], HSFs[z-1] = LSF, HSF
            HSFs_rms.append(get_rms(HSF, verbose=verbose))
            print(z, end=', ')
        spsd = [rms**2 for rms in HSFs_rms]
        spsd = [spsd[0]] + spsd
        fits.writeto(spsd_name%'spsd', np.float32(spsd))
        fits.writeto(spsd_name%'LSFs', np.float32(LSFs))
        fits.writeto(spsd_name%'HSFs', np.float32(HSFs))
    return spsd

def psd_temporal(zpols, pol_order=3, win_size=51, savgol=True):
    '''
    pol_order: filter polynomial order
    win_size: filter window size
    '''
    N = len(zpols)      # time domain
    M = N // 2 + 1      # frequency domain
    yf = np.abs(np.array([np.fft.fft(yt)[:M] for yt in zpols.swapaxes(0,1)]))
    if savgol is True:  # savgol filter 
        warnings.filterwarnings("ignore")
        yf = savgol_filter(yf, win_size, pol_order)
    tpsd = yf**2
    return tpsd