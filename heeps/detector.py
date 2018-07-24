import proper
from astropy.io import fits 

def detector(wfo,f_lens,nd,coronagraph_type,prefix,Debug=False):
    n = proper.prop_get_gridsize(wfo)
    if (n >= nd):
        proper.prop_propagate(wfo, f_lens, "to reimaging lens")
        proper.prop_lens(wfo, f_lens, "apply reimaging lens")
        proper.prop_propagate(wfo, f_lens, "final focus")
        (wfo, sampling) = proper.prop_end(wfo, NOABS = False) # conclude the simulation --> noabs= the wavefront array will be complex
    else: 
        print('Error: final image is bigger than initial grid size')            
    psf = wfo[int(n/2-nd/2):int(n/2+nd/2),int(n/2-nd/2):int(n/2+nd/2)]
    out_dir = str('./output_files/')
    if (Debug==True): 
        fits.writeto(out_dir + prefix + '_' + coronagraph_type+'_PSF'+'.fits', psf, overwrite=True)
    return psf 

