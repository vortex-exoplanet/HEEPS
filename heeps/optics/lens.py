import proper

def lens(wf, focal=660, offset_before=0, offset_after=0, offset_light_trap=0, 
        diam_light_trap=0, **conf):

    # propagation before lens
    proper.prop_propagate(wf, focal + offset_before)
    # Fourier transform of an image using a lens
    proper.prop_lens(wf, focal)

    # check for light trap
    if offset_light_trap == 0 and diam_light_trap == 0:
        # propagation after lens
        proper.prop_propagate(wf, focal + offset_after)
    else:
        # add light trap
        proper.prop_propagate(wf, focal - offset_light_trap)
        proper.prop_circular_aperture(wf, diam_light_trap, 0, 0, NORM=True)
        proper.prop_propagate(wf, offset_light_trap + offset_after)