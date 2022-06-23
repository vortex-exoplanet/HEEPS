from .load_errors import load_errors
from .propagate_cube import propagate_cube

def propagate(wf, nframes=10, nstep=1, nframes_avg=10, avg=False, tag=None,
        onaxis=True, send_to=None, savefits=False, verbose=False, **conf):

    if verbose is True:
        print('%s-axis PSF:'%{True:'On',False:'Off'}[onaxis])
        print('\u203e'*{True:12,False:13}[onaxis])

    if avg is True:
        conf.update(
            nframes = nframes_avg,
            nstep = 1,
            send_to = None,
        )
    else:
        conf.update(nframes=nframes, nstep=nstep, send_to=send_to)

    phase_screens, amp_screens, tiptilts, apo_misaligns, ls_misaligns = \
        load_errors(verbose=verbose, **conf)

    psfs = propagate_cube(wf, phase_screens=phase_screens,
        amp_screens=amp_screens, tiptilts=tiptilts, apo_misaligns=apo_misaligns,
        ls_misaligns=ls_misaligns, tag=tag, onaxis=onaxis, avg=avg, 
        savefits=savefits, verbose=verbose, **conf)

    return psfs