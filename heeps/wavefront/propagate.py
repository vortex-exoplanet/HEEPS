from .load_errors import load_errors
from .propagate_cube import propagate_cube

def propagate(wf, tag=None, onaxis=True, avg=False, send_to=None, savefits=False, verbose=False, **conf):

    if verbose is True:
        print('%s-axis PSF:'%{True:'On',False:'Off'}[onaxis])
        print('\u203e'*{True:12,False:13}[onaxis])

    phase_screens, amp_screens, tiptilts, apo_misaligns, ls_misaligns = \
        load_errors(verbose=verbose, **conf)

    psfs = propagate_cube(wf, phase_screens=phase_screens,
        amp_screens=amp_screens, tiptilts=tiptilts, apo_misaligns=apo_misaligns,
        ls_misaligns=ls_misaligns, tag=tag, onaxis=onaxis, avg=avg, 
        send_to=send_to, savefits=savefits, verbose=verbose, **conf)

    return psfs