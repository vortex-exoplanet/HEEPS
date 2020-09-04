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
    onaxis = True,
    cpu_count = None
)
conf = heeps.config.read_config(verbose=True, **conf)

for conf['band'] in conf['bands']:
    for conf['mode'] in conf['modes']:
    
        # 2. (optional) Update config parameters. ATTENTION, the following parameters 
        # will be updated to match the selected spectral band and HCI mode:
        #   lam, pscale, flux_star, flux_bckg, ls_dRspi, ls_dRint, npupil, ndet,
        #   ravc_t, ravc_r
        conf = heeps.config.update_config(verbose=True, **conf) 

        # 3. Load wavefront errors
        phase_screens, pointing_errs, apo_drifts = heeps.wavefront.load_errors(verbose=True, **conf)

        # 4. Load entrance pupil, and create 'wavefront' object
        wf, pup = heeps.pupil.pupil(verbose=True, **conf)

        # 5. Propagate cube of psfs
        psfs = heeps.wavefront.propagate_cube(wf, conf, phase_screens=phase_screens, \
            savefits=True)

print('Simulation finished.')

#########

if False:
    import matplotlib.pyplot as plt
    plt.imshow(pup)
    