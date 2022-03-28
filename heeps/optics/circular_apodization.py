import proper

def circular_apodization(wf, radius, t_in, t_out, xc=0, yc=0, NORM=False):

    apodizer = proper.prop_ellipse(wf, radius, radius, xc, yc, NORM=NORM)

    if (t_in > t_out):
        apodizer = apodizer*(t_in - t_out) + t_out
    else:
        apodizer = apodizer*(t_out - t_in) + t_in

    return apodizer
