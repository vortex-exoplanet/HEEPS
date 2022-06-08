import numpy as np

def paralang(npts, dec, lat, duration=3600):
    """ parallactic angle """
    # hour angle to deg conversion
    ha = duration/3600*360/24
    # angles in rad
    hr = np.deg2rad(np.linspace(-ha/2, ha/2, npts))
    dr = np.deg2rad(dec)
    lr = np.deg2rad(lat)
    # parallactic angle in deg
    pa = -np.rad2deg(np.arctan2(-np.sin(hr), 
                                 np.cos(dr)*np.tan(lr) - np.sin(dr)*np.cos(hr)))
    return pa