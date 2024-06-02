import heeps
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = True


def make_heeps_app_psf(ngs_mag=3.5, filter_name="L", dit=0.3,
                       app_phase_ramp_params={"offset": 10, "angle": 30}):

    # set up the HEEPS config
    conf = heeps.config.read_config(verbose=False)      # just a dict
    conf.update(mode="APP", mag=ngs_mag, band=filter_name, dit=dit,
                app_phase_ramp_params=app_phase_ramp_params,
                nframes=4, disp=False, ndet=1024, ngrid=1024)

    # Make a pupil wavefront array
    wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)

    # Generate the 3 APP PSF components
    leakage_psf = heeps.wavefront.propagate(wf, onaxis=False, avg=True, **conf)
    left_half_psf = heeps.wavefront.propagate(wf, onaxis=True, avg=True,**conf)
    conf["mode"] = "APPneg"
    right_half_psf = heeps.wavefront.propagate(wf, onaxis=True, avg=True, **conf)

    # Combine the 3 components
    psfs = np.stack([left_half_psf, leakage_psf, right_half_psf], axis=0)
    psf_lobe_weight = conf["app_single_psf"]
    psf_weights = [psf_lobe_weight, 1 - 2 * psf_lobe_weight, psf_lobe_weight]
    psfs *= np.array(psf_weights)[:, None, None]
    psf = psfs.sum(axis=0)
    psf /= psf.sum()            # Originally this was psf.max() @gotten, why?
    assert abs(psf.sum() - 1) < 0.02

    return psf


psf = make_heeps_app_psf()

if PLOTS:
    plt.imshow(psf, norm=LogNorm(vmin=1e-5), origin="lower")
    plt.colorbar()
    plt.show()

# fits.writeto("off_axis_psf.fits", off_axis_psf, overwrite=True)
