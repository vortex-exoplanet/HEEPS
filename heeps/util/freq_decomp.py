import heeps.util.img_processing as impro
import astropy.convolution as astroconv
import numpy as np
import warnings

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
    kernel = impro.resize_img(kernel, nkernel)
    kernel /= np.sum(kernel)                        # need to normalize the kernel
    
    return kernel

def spatial(allSF, kernel, npupil=None, norm=False, verbose=False):
    
    # mask with nans
    allSF[allSF==0] = np.nan
    # get low and high spatial frequencies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # NANs
        LSF = astroconv.convolve(allSF, kernel, boundary='extend')
        HSF = allSF - LSF
    # print rms
    if verbose is True:
        print('rms(all SF) = %3.2f'%(np.nanstd(allSF)))
        print('rms(LSF) = %3.2f'%(np.nanstd(LSF)))
        print('rms(HSF) = %3.2f'%(np.nanstd(HSF)))
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
        allSF = impro.resize_img(allSF, npupil)
        LSF = impro.resize_img(LSF, npupil)
        HSF = impro.resize_img(HSF, npupil)
    
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

