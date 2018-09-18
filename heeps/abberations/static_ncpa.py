import numpy as np
import cv2
import proper

def static_ncpa(wf, npupil, screen):
    
    n = int(proper.prop_get_gridsize(wf))

    screen_pixels = (screen.shape)[0] # size of the phase screen
    sf = float(npupil)/float(screen_pixels) # scaling factor the phase screen to the simulation pupil size

    screen_scale = cv2.resize(screen.astype(np.float32), (0,0), fx = sf, fy = sf, interpolation=cv2.INTER_LINEAR) # scale the the phase screen to the simulation pupil size
    screen_large = np.zeros((n,n)) # define an array of n-0s, where to insert the screen
    screen_large[int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2),int(n/2)+1-int(npupil/2)-1:int(n/2)+1+int(npupil/2)] =screen_scale # insert the scaled screen into the 0s grid

    proper.prop_add_phase(wf, screen_large)
    
    return wf
