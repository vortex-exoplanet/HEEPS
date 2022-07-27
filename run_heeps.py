#!/usr/bin/env python3

# ===================== #
# End-to-end simulation #
# ===================== #

import heeps
import sys

def run_heeps(savename='contrast_curve.png'):

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

    # 6. Produce a raw contrast curve
    sep, raw = heeps.contrast.cc_raw(savefits=True, verbose=True, **conf)

    # 7. Produce a 5-sigma sensitivity (contrast) curve
    sep1, adi1 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)

    # 8. Add star flux + background flux + photon noise
    conf['add_bckg'] = True
    sep2, adi2 = heeps.contrast.cc_adi(savepsf=True, savefits=True, verbose=True, **conf)

    # 9. Create a figure 
    if savename != '':
        xlabel = 'Angular separation $[arcsec]$'
        ylabel_adi = '5-$\sigma$ sensitivity (contrast)'
        ylabel_raw = 'raw contrast'
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 6))
        fig.subplots(2, 1, sharex=True)
        fig.subplots_adjust(hspace=0.02)
        axes = fig.axes
        axes[0].set_ylim(1e-7, 1e-2)
        axes[1].set_ylim(1e-7, 9e-3)
        axes[1].set_xlabel(xlabel)
        for j, (ax, ylabel) in enumerate(zip(axes, [ylabel_raw, ylabel_adi])):
            ax.set_ylabel(ylabel)
            ax.grid(True), ax.grid(which='minor', linestyle=':')
            ax.loglog()
            ax.xaxis.set_major_formatter(plt.ScalarFormatter())
            ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
            ax.set_xlim(0.02, 0.75)
        axes[0].plot(sep, raw, 'C0', label='RAVC', marker='d', markevery=0.12, markersize=4)
        axes[0].legend()
        axes[0].set_title('HCI modes; scao only; star mag L = 6')
        axes[1].plot(sep1, adi1, 'C0', label='RAVC', marker='d', markevery=0.12, markersize=4)
        axes[1].plot(sep2, adi2, ':C0', label='background', marker='d', markevery=0.12, markersize=4)
        axes[1].legend(ncol=2, loc='upper right')
        fig.savefig('%s/%s'%(conf['dir_output'], savename), dpi=300, transparent=True);

if __name__ == "__main__":
    '''
    Terminal command line example
    > python run_heeps.py test.png
    '''
    if len(sys.argv) > 1:
        run_heeps(sys.argv[1])
    else:
        run_heeps()