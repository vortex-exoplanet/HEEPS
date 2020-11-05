import matplotlib.pyplot as plt
import proper
import numpy as np

def imshow_amp(wf, npupil=None, margin=0):

    amp = proper.prop_get_amplitude(wf)
    ngrid = amp.shape[0]
    if npupil is None:
        npupil = ngrid
    start = int((ngrid - npupil + npupil%2)/2 - margin)
    start = 0 if start < 0 else start
    end = ngrid - start + npupil%2
    plt.imshow(amp[start:end,start:end], origin='lower')

def imshow_phi(wf, npupil=None, margin=0):

    amp = proper.prop_get_phase(wf)
    ngrid = amp.shape[0]
    if npupil is None:
        npupil = ngrid
    start = int((ngrid - npupil + npupil%2)/2 - margin)
    start = 0 if start < 0 else start
    end = ngrid - start + npupil%2
    plt.imshow(amp[start:end,start:end], origin='lower')