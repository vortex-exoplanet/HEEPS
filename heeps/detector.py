import proper

def detector(wfo,f_lens,nd):
    n = proper.prop_get_gridsize(wfo)
    if (n >= nd):
        proper.prop_propagate(wfo, f_lens, "to reimaging lens")
        proper.prop_lens(wfo, f_lens, "apply reimaging lens")
        proper.prop_propagate(wfo, f_lens, "final focus")
        (wfo, sampling) = proper.prop_end(wfo, NOABS = False) # conclude the simulation --> noabs= the wavefront array will be complex
    else: 
        print('Error: final image is bigger than initial grid size')            
    return wfo[int(n/2-nd/2):int(n/2+nd/2),int(n/2-nd/2):int(n/2+nd/2)]

