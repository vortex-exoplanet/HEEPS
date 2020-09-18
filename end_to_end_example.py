#!/usr/bin/env python3

# ============================= #
# End-to-end simulation example #
# ============================= #

import heeps

# 1. Create a config dictionary with your simulation parameters. 
#    Undefined parameters get default values from calling read_config
conf = dict(
    dir_current = '$HOME/Desktop',
    bands = ['L'],
    modes = ['RAVC','CVC'],
    cpu_count = None,
    add_scao=True
)
conf = heeps.config.read_config(**conf)

for conf['band'] in conf['bands']:
    for conf['mode'] in conf['modes']:
    
        # 2. Update config parameters. The following parameters 
        # will be updated to match the selected spectral band and HCI mode:
        #   lam, pscale, flux_star, flux_bckg, ls_dRspi, ls_dRint, npupil, ndet,
        #   ravc_t, ravc_r
        conf = heeps.config.update_config(verbose=True, **conf) 

        # 3. Load wavefront errors
        phase_screens, pointing_errs, apo_drifts = heeps.wavefront.load_errors(verbose=True, **conf)

        # 4. Load entrance pupil, and create 'wavefront' object
        wf = heeps.pupil.pupil(verbose=True, savefits=True, **conf)

        # 5. Propagate one frame of offaxis psf (i.e. planet)
        heeps.wavefront.propagate_one(wf, savefits=True, onaxis=False, **conf)

        # 6. Propagate cube of onaxis psfs (i.e. star)
        heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
            verbose=True, savefits=True, **conf)

print('Simulation finished.')
