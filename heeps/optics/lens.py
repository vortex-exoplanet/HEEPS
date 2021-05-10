import proper

def lens(wf, focal=660, offset_before=0, offset_after=0, **conf):
    
    # propagation before lens
    proper.prop_propagate(wf, focal + offset_before)
    # Fourier transform of an image using a lens
    proper.prop_lens(wf, focal)
    # propagation after lens
    proper.prop_propagate(wf, focal + offset_after)