import proper

def lens(wf, focal=660, offset_before=0, offset_after=0, offset_light_trap=0, 
        diam_light_trap=0, **conf):

    # check for light trap
    if diam_light_trap == 1:
        print('no light trap')
        proper.prop_propagate(wf, offset_light_trap)
    elif offset_light_trap != 0:
        print('light trap')
        proper.prop_propagate(wf, offset_light_trap)
        proper.prop_circular_aperture(wf, diam_light_trap, 0, 0, NORM=True)
        proper.prop_propagate(wf, focal - offset_light_trap)
        proper.prop_lens(wf, focal)
        proper.prop_propagate(wf, focal)
    else:
        # propagation before lens
        proper.prop_propagate(wf, focal + offset_before)
        # Fourier transform of an image using a lens
        proper.prop_lens(wf, focal)
        # propagation after lens
        proper.prop_propagate(wf, focal + offset_after)