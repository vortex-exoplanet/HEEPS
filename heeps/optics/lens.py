import proper
import numpy as np

def lens(wf, focal):
    proper.prop_propagate(wf, focal)
    proper.prop_lens(wf, focal)
    proper.prop_propagate(wf, focal)

def lens2(wf, focal):
    proper.prop_propagate(wf, focal)
    wf._wfarr = np.fft.fft2(wf._wfarr)/wf._ngrid
    proper.prop_propagate(wf, focal)

def lens_offset(wf, focal, offset_before=0, offset_after=0):
    proper.prop_propagate(wf, focal + offset_before)
    proper.prop_lens(wf, focal)
    proper.prop_propagate(wf, focal + offset_after)