import numpy as np

def paralang(npts, dec, lat, duration=3600):
    """ parallactic angle 
    
    Examples:
        np.sum(abs(paralang(2, -61, -24.59, 3600)))
        22.870920714712355

        np.sum(abs(paralang(2, -61, -24.59, 3600*5)))
        104.42140816101042
        
    Returns:
        pa (ndarray): parallactic angles in deg
    """

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

def field_rotation(pa):
    """ 
    total field rotation in deg 
    
    PARAMETERS:
    pa (ndarray): parallactic angles in deg
    """
    return np.ptp(np.unwrap(pa, period=360))