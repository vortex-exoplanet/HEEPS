import proper
import numpy.fft as fft
import pyfftw.interfaces.numpy_fft as fftw

def lens(wf, focal=660, offset_before=0, offset_after=0, **conf):
    
    # propagation before lens
    proper.prop_propagate(wf, focal + offset_before)
    # Fourier transform of an image using a lens
    proper.prop_lens(wf, focal)
    # propagation after lens
    proper.prop_propagate(wf, focal + offset_after)

def lens_test(wf, focal=660, lens_method='proper', offset_before=0, offset_after=0, **conf):
    
    # propagation before lens
    proper.prop_propagate(wf, focal + offset_before)
    # Fourier transform of an image using a lens
    if lens_method == 'proper':
        proper.prop_lens(wf, focal)
    elif lens_method == 'numpy':
        wf._wfarr = fft.fft2(wf._wfarr)/wf._ngrid
    elif lens_method == 'pyfftw':
        wf._wfarr = fftw.fft2(wf._wfarr)/wf._ngrid
    # propagation after lens
    proper.prop_propagate(wf, focal + offset_after)