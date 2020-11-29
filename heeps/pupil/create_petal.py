from heeps.util.polar_coord import ri_ti
import numpy as np

def create_petal(select_petal, npetals, npupil):

    # petal start and end angles
    pet_angle = 2*np.pi/npetals
    pet_start = pet_angle/2 + (select_petal - 1)*pet_angle
    pet_end = pet_start + pet_angle
    # grid angles
    ri, ti = ri_ti(npupil)
    # petal angles must be 0-2pi
    ti %= (2*np.pi)
    pet_start %= (2*np.pi)
    pet_end %= (2*np.pi)
    # check if petal crosses 0 angle
    if pet_end - pet_start < 0:
        pet_start = (pet_start + np.pi)%(2*np.pi) - np.pi
        pet_end = (pet_end + np.pi)%(2*np.pi) - np.pi
        ti = (ti + np.pi)%(2*np.pi) - np.pi
    # create petal and add to grid
    petal = np.uint8((ti>=pet_start) * (ti<=pet_end))
    
    return petal