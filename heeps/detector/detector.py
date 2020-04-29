import proper
from astropy.io import fits 

def detector(wfo, conf):
    f_lens = conf['focal']
    ndet = conf['ndet']
    gridsize = conf['gridsize']
    mode = conf['mode']
    prefix = conf['prefix']
    
    if (gridsize >= ndet):
        proper.prop_propagate(wfo, f_lens, "to reimaging lens")
        proper.prop_lens(wfo, f_lens, "apply reimaging lens")
        proper.prop_propagate(wfo, f_lens, "final focus")
        # conclude the simulation --> noabs = wavefront array will be complex
        (wfo, sampling) = proper.prop_end(wfo, NOABS = False)
    else: 
        print('Error: final image is bigger than initial grid size')
    start = int(gridsize/2 - ndet/2) + 1
    end = int(gridsize/2 + ndet/2)
    psf = wfo[start:end, start:end]
    out_dir = str('./output_files/')
    
    return psf
