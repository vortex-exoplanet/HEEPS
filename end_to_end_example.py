#!/usr/bin/env python3

# ============================= #
# End-to-end simulation example #
# ============================= #

import heeps

# 1. Create a config dictionary with your simulation parameters. 
#    Undefined parameters get default values from calling read_config
conf = dict(
    bands = ['L'],
    modes = ['RAVC','CVC'],
    dir_current = '$HOME/Desktop/heeps_analysis',
    dir_output = 'scao_only',
    add_scao = True,
    send_to = None, #'cdelacroix@uliege.be',
    cpu_count = None
)
conf = heeps.config.read_config(verbose=False, **conf)

for conf['band'] in conf['bands']:
    for conf['mode'] in conf['modes']:
    
        # 2. Update config parameters. The following parameters 
        # will be updated to match the selected spectral band and HCI mode:
        #   lam, pscale, flux_star, flux_bckg, ls_dRspi, ls_dRint, npupil, ndet,
        #   ravc_t, ravc_r
        conf = heeps.config.update_config(saveconf=True, verbose=True, **conf) 

        # 3. Load wavefront errors
        phase_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)

        # 4. Load entrance pupil, and create 'wavefront' object
        wf = heeps.pupil.pupil(savefits=True, verbose=True, **conf)

        # 5. Propagate one frame of offaxis psf (i.e. planet)
        psf = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=False, **conf)

        # 6. Propagate cube of onaxis psfs (i.e. star)
        psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
            tiptilts=tiptilts, misaligns=misaligns, onaxis=True, savefits=True, \
            verbose=True, **conf)

# 7. Produce 5-sigma sensitivity (contrast) curves for each set of PSFs (modes, bands)
if True:
    sep, sen = heeps.contrast.adi_one(savepsf=True, savefits=True, verbose=True, **conf)
