import proper
from astropy.io import fits

def lens(wf, focal=660, zoffset_before=0, zoffset_after=0,
        dist_start=0, diam_start=0, dist_end=0, diam_end=0,
        lt_write_name=None, **conf):

    # check for light trap BEFFORE lens
    if dist_start == 0 and diam_start == 0:
        # propagation before lens
        proper.prop_propagate(wf, focal + zoffset_after)
    else:
        # add light trap
        proper.prop_propagate(wf, dist_start)
        if type(diam_start) is tuple:
            proper.prop_rectangular_aperture(wf, diam_start[0], diam_start[1], 0, 0, NORM=True)
        else:
            proper.prop_circular_aperture(wf, diam_start, 0, 0, NORM=True)
        if lt_write_name is not None:
            fits.writeto(lt_write_name, proper.prop_get_amplitude(wf), overwrite=True)
        proper.prop_propagate(wf, focal - dist_start + zoffset_after)
    
    # Fourier transform of an image using a lens
    proper.prop_lens(wf, focal)

    # check for light trap AFTER lens
    if dist_end == 0 and diam_end == 0:
        # propagation after lens
        proper.prop_propagate(wf, focal + zoffset_before)
    else:
        # add light trap
        proper.prop_propagate(wf, focal - dist_end)
        if type(diam_end) is tuple:
            proper.prop_rectangular_aperture(wf, diam_end[0], diam_end[1], 0, 0, NORM=True)
        else:
            proper.prop_circular_aperture(wf, diam_end, 0, 0, NORM=True)
        if lt_write_name is not None:
            fits.writeto(lt_write_name, proper.prop_get_amplitude(wf), overwrite=True)
        proper.prop_propagate(wf, dist_end + zoffset_before)