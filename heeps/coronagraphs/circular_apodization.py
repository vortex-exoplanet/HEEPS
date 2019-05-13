import proper

def circular_apodization(wf, radius, t_in, t_out, xc = 0.0, yc = 0.0, **kwargs):

    if ("NORM" in kwargs and kwargs["NORM"]):
        norm = True
    else:
        norm = False

    if (t_in > t_out):
        apodizer = proper.prop_shift_center(proper.prop_ellipse(wf, radius, \
                    radius, xc, yc, NORM = norm))*(t_in-t_out)+t_out
    else:
        apodizer = proper.prop_shift_center(proper.prop_ellipse(wf, radius, \
                    radius, xc, yc, NORM = norm))*(t_out-t_in)+t_in

    return apodizer
