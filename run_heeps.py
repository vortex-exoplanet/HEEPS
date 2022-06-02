#!/usr/bin/env python3

# ===================== #
# End-to-end simulation #
# ===================== #

import heeps
import sys

def run_heeps(savename='ADI_contrast_curve.png'):

    # 1. Create a config dictionary with your simulation parameters in read_config
    conf = heeps.config.read_config(verbose=False)

    # 2. Update config parameters. The following parameters 
    # will be updated to match the selected spectral band and HCI mode:
    #   lam, pscale, flux_star, flux_bckg, ls_dRspi, ls_dRint, npupil, ndet,
    #   ravc_t, ravc_r
    conf = heeps.config.update_config(saveconf=True, verbose=True, **conf) 

    # 3. Load entrance pupil, and create 'wavefront' object
    wf = heeps.pupil.pupil(savefits=True, verbose=True, **conf)

    # 4.  Create a PSF template for injecting fake exoplanets
    heeps.wavefront.propagate(wf, onaxis=False, avg=True, savefits=True, verbose=True, **conf);

    # 5. Create a cube of on-axis psfs (=star)
    heeps.wavefront.propagate(wf, onaxis=True, savefits=True, verbose=True, **conf);

    # 6. Produce a 5-sigma sensitivity (contrast) curve
    sep1, sen1 = heeps.contrast.adi_one(savepsf=True, savefits=True, verbose=True, **conf)

    # 7. Add star flux + background flux + photon noise
    conf['add_bckg'] = True
    sep2, sen2 = heeps.contrast.adi_one(savepsf=True, savefits=True, verbose=True, **conf)

    # 8. Create a figure 
    if savename != '':
        import matplotlib.pyplot as plt
        from matplotlib.pylab import ScalarFormatter
        plt.figure(figsize=(9.5, 4))
        plt.grid(True), plt.grid(which='minor', linestyle=':')
        plt.loglog(), plt.gca().xaxis.set_major_formatter(ScalarFormatter())
        plt.xlabel('Angular separation $[arcsec]$')
        plt.ylabel('5-$\sigma$ sensitivity (contrast)')
        plt.title('%s-band %s, SCAO only, %s frames'%(conf['band'], conf['mode'], conf['nframes']))
        plt.plot(sep1, sen1, 'C0', label='w/o background', marker='d', markevery=0.12, markersize=4)
        plt.plot(sep2, sen2, 'C0:', label='with background (star mag=5)', marker='d', markevery=0.12, markersize=4)
        plt.legend()
        plt.xlim(0.02, 0.75)
        plt.ylim(1e-8, 1e-2)
        plt.xticks([0.02, 0.05, 0.1, 0.2, 0.5])
        plt.savefig('%s/%s'%(conf['dir_output'], savename), dpi=300, transparent=True);

if __name__ == "__main__":
    '''
    Terminal command line example
    > python run_heeps.py test.png
    '''
    if len(sys.argv) > 1:
        run_heeps(sys.argv[1])
    else:
        run_heeps()