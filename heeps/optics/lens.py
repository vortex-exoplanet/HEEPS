import proper

def lens(wf, focal=660, vc_zoffset_before=0, vc_zoffset_after=0,
        lt_dist_start=0, lt_diam_start=0, lt_dist_end=0, lt_diam_end=0, **conf):

    # check for light trap BEFFORE lens
    if lt_dist_start == 0 and lt_diam_start == 0:
        # propagation before lens
        proper.prop_propagate(wf, focal + vc_zoffset_after)
    else:
        # add light trap
        proper.prop_propagate(wf, lt_dist_start)
        if type(lt_diam_start) is tuple:
            proper.prop_rectangular_aperture(wf, lt_diam_start[0], lt_diam_start[1], 0, 0, NORM=True)
        else:
            proper.prop_circular_aperture(wf, lt_diam_start, 0, 0, NORM=True)
        proper.prop_propagate(wf, focal - lt_dist_start + vc_zoffset_after)
    
    # Fourier transform of an image using a lens
    proper.prop_lens(wf, focal)

    # check for light trap AFTER lens
    if lt_dist_end == 0 and lt_diam_end == 0:
        # propagation after lens
        proper.prop_propagate(wf, focal + vc_zoffset_before)
    else:
        # add light trap
        proper.prop_propagate(wf, focal - lt_dist_end)
        if type(lt_diam_end) is tuple:
            proper.prop_rectangular_aperture(wf, lt_diam_end[0], lt_diam_end[1], 0, 0, NORM=True)
        else:
            proper.prop_circular_aperture(wf, lt_diam_end, 0, 0, NORM=True)
        proper.prop_propagate(wf, lt_dist_end + vc_zoffset_before)