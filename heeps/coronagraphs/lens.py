import proper

def lens(wfo, f_lens):
    proper.prop_propagate(wfo, f_lens)
    proper.prop_lens(wfo, f_lens)
    proper.prop_propagate(wfo, f_lens)
