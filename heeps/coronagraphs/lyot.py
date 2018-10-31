import numpy as np
import proper
import os
from ..fits import writefield 
from ..fits import readfield





def lyot(wfo, conf):
    tmp_dir = conf['TMP_DIR']
    n = int(proper.prop_get_gridsize(wfo))
    ofst = 0 # no offset
    
    f_lens = conf['F_LENS']
    diam = conf['DIAM']
    pixelsize = conf['PIXEL_SCALE']
    Debug_print = conf['DEBUG_PRINT'] 
    diam_CL = conf['DIAM_CL'] # classical lyot diam in lam/D
    
    wavelength = proper.prop_get_wavelength(wfo) 
    gridsize = proper.prop_get_gridsize(wfo)
    beam_ratio = pixelsize*4.85e-9/(wavelength/diam)
    calib = str(int(beam_ratio*100))+str('_')+str(gridsize)
    my_file = str(tmp_dir+'zz_perf_'+calib+'_r.fits')

    proper.prop_propagate(wfo, f_lens, 'inizio') # propagate wavefront
    proper.prop_lens(wfo, f_lens, 'focusing lens LC') # apply lens
    proper.prop_propagate(wfo, f_lens, 'LC') # propagate wavefront

    if (os.path.isfile(my_file)==True):
        vvc = readfield(tmp_dir,'zz_vvc_'+calib) # read the theoretical vortex field
        vvc = proper.prop_shift_center(vvc)
        scale_psf = wfo._wfarr[0,0]
        psf_num = readfield(tmp_dir,'zz_psf_'+calib) # read the pre-vortex field
        psf0 = psf_num[0,0]
        psf_num = psf_num/psf0*scale_psf
        perf_num = readfield(tmp_dir,'zz_perf_'+calib) # read the perfect-result vortex field
        perf_num = perf_num/psf0*scale_psf
        wfo._wfarr = (wfo._wfarr - psf_num)*vvc + perf_num # the wavefront takes into account the real pupil with the perfect-result vortex field
    else:
        samp = 100
        n_samp = samp*int(np.ceil(diam_CL/beam_ratio))
        Rext = round(samp/2.*diam_CL/beam_ratio)
        center = (n_samp - 1)/2. # 513
        temp = np.zeros([n_samp,n_samp])
        for i in range(n_samp):
            for j in range(n_samp):
                r = np.sqrt((i - center)**2 + (j - center)**2)
                if r > Rext:
                    temp[i,j] = 1
        
        import scipy.misc
        n_samp2 = int(n_samp/samp)
        mask_samp = scipy.misc.imresize(temp, (n_samp2, n_samp2), interp="bilinear")/255.
        
        npts = 1024
        mask = np.ones([npts,npts])
        mask[int(np.ceil((npts - n_samp2)/2.)):int(np.ceil((npts + n_samp2)/2.)), \
                int(np.ceil((npts - n_samp2)/2.)):int(np.ceil((npts + n_samp2)/2.))]\
                = mask_samp
        
        mask = proper.prop_shift_center(mask)
        wfo._wfarr *= mask # apply lyot mmask
    
    proper.prop_propagate(wfo, f_lens, "propagate to pupil reimaging lens")  
    proper.prop_lens(wfo, f_lens, "apply pupil reimaging lens")
    proper.prop_propagate(wfo, f_lens, "lyot stop")
            
    return wfo






