#!/usr/bin/env python3

# ============================= #
# End-to-end simulation example #
# ============================= #

import heeps

# 1. Create a config dictionary with your simulation parameters. 
#    Undefined parameters get default values from calling read_config
conf = dict(
    band = 'L',
    mode = 'RAVC',
    cpu_count = 10,
)
conf = heeps.config.read_config(verbose=True, **conf)

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
psf = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)

# 6. Propagate cube of onaxis psfs (i.e. star)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    tiptilts=tiptilts, misaligns=misaligns, onaxis=True, savefits=True, \
    verbose=True, **conf)

# 7. Produce 5-sigma sensitivity (contrast) curves for each set of PSFs (modes, bands)
sep, sen = heeps.contrast.adi_one(savepsf=True, savefits=True, verbose=True, **conf)

# 8. Create a figure 
if True:
    %pylab
    %matplotlib inline
    figure(figsize(12,4))
    grid(True), grid(which='minor', linestyle=':')
    loglog(), gca().xaxis.set_major_formatter(ScalarFormatter())
    xlabel('Angular separation $[arcsec]$')
    ylabel('5-$\sigma$ sensitivity (contrast)')
    label = '%s-band %s'%(conf['band'], conf['mode'])
    plot(sep, sen, label=label, marker='d', markevery=0.12, markersize=4)
    legend()
    xticks([0.02, 0.1, 0.5, 1])
    xlim(0.02,1)
    ylim(1e-6,1e-2)
    savefig('example_adi.png', dpi=300, transparent=True)