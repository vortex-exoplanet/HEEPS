"""
HEEPS disc injection utilities.

This module provides functions to inject a synthetic disc model into HEEPS
PSF cubes for later post‑processing with VIP.

The entry point is :func:`create_disc_sequence`, which creates a mock
observation sequence with a user-provided disc model.

"""

import numpy as np
from typing import Union
import os

from vip_hci.fits import open_fits
from vip_hci.preproc import frame_crop, cube_crop_frames, frame_shift
from vip_hci.fm import cube_inject_fakedisk

from heeps.util.psf_template import psf_template
from heeps.util.paralang import paralang

__all__ = ['create_disc_sequence']
__author__ = "Iain Hammond"

def create_disc_sequence(
        disc_model : np.ndarray,
        db_path : str = "vortex_psf_grid/",
        seeing_q : int = 2,
        source_xy : Union[tuple, list, None] = None,
        extinction : Union[int, float] = 0,
        starphot : Union[float, int, None] = 1e11,
        do_rdi : bool = True,
        rdi_mag : Union[int, float, None] = None,
        rdi_duration : Union[int, float, None] = None,
        imlib : str = "opencv",
        **conf
):
    """
    Create a mock observation sequence by injecting a synthetic disc model into a HEEPS PSF. The function will attempt
    to load from the PSF grid directory based on the conf parameters.

    The code assumes that the disc model has the star included to convert to units of contrast.

    Parameters
    ----------
    disc_model : np.ndarray
        2‑D array representing the disc model to be injected. NaNs are replaced with zeros. Odd-sized models with the
        stellar flux on the centre pixel are preferred.
    db_path : str
        Path to the PSF grid directory (the folder that contains run subdirectories).
    seeing_q : int, default=2
        Seeing quartile (1-3) to select the PSF. Uses a default seeing of Q2 if not specified.
    source_xy : Union[tuple, list, None], optional
        (x, y) pixel coordinates of the source for which extinction is applied.
    extinction : Union[int, float], default=0
        Extinction in magnitudes to apply to the source at ``source_xy``.
    starphot : float, default=1e11
        Photometric scaling factor for the star.
    do_rdi : bool, default=True
        Whether to prepare a reference cube for RDI (Reference Differential Imaging). If ``True``, a reference cube is
        returned as the third output of this function.
    rdi_mag : Union[int, float, None], optional
        Magnitude to use for the RDI reference. If ``None``, magnitude in conf is used (i.e., the same as the science
        target).
    rdi_duration : Union[int, float, None], optional
        Desired duration (in seconds) of the RDI reference sequence. If ``None``, a default of 20% of the on‑axis
        sequence length is used.
    imlib : str, default="opencv"
        Image library to use for image processing. Use "opencv" to go fast, or "vip-fft" for better flux conservation
        at the cost of speed. See VIP documentation for details.
    conf : dict, optional
        Configuration dictionary with parameters from HEEPS (e.g., band, dit, mag, etc.).

    Returns
    -------
    psf_ON : numpy.ndarray
        Sequence with the fake disc injected into the on‑axis PSF cube.
    pa : numpy.ndarray
        Array of parallactic angles used for the injection.
    psf_RDI : numpy.ndarray, optional
        Reference cube for RDI, if ``do_rdi`` is ``True``. Otherwise, this output is not returned.
    """
    # if mag is not in the grid, round to the closest available magnitude
    conf["mag"] = _check_grid(conf["mag"], conf["band"], conf["mode"])

    if do_rdi and rdi_mag is None:
        rdi_mag = conf["mag"]
        print("Note: rdi_mag was not set. Using mag in the conf dictionary.", flush=True)

    # prepare the path for loading the PSF grid directory
    # make sure db_path ends with a slash
    if not db_path.endswith("/"):
        db_path += "/"

    # new grid folder structure and naming convention for the April 2026 grid files,
    # for example, vortex_psf_grid/L_CVC_s=Q2_mag=8.5_3600s_100ms
    grid_dir = os.path.join(
        db_path,
        f"{conf['band']}_{conf['mode']}_s=Q{seeing_q}_mag={conf['mag']}_{conf['duration']}s_{int(conf['dit'] * 1000)}ms"
    )
    loadname = os.path.join(grid_dir, '%s_PSF_bckg1_%s_%s.fits' % ('%s', conf['band'], conf['mode']))

    # set nans to zero
    disc_model = np.nan_to_num(disc_model, nan=0)

    # remove any extra dimensions (e.g., MCFOST cubes)
    # first, squeeze singleton dimensions
    disc_model = np.squeeze(disc_model)
    # then collapse any remaining leading dimensions to a 2-D image
    while disc_model.ndim > 2:
        disc_model = disc_model[0]

    # record the stellar flux in the model for flux normalisation, and set its pixel to 0 (it will be in the on-axis PSF)
    # some codes have the star on the centre 4 pixels instead of the centre pixel only so we need to sum the values
    mask = disc_model == disc_model.max()
    star_val = disc_model[mask].sum()
    disc_model[mask] = 0

    # ensure the model has an odd size
    if disc_model.shape[-1] % 2 == 0:
        print("Converting input disc model to odd-dimensions.", flush=True)
        disc_model = frame_shift(disc_model, shift_y=0.5, shift_x=0.5, imlib=imlib)
        disc_model = disc_model[1:, 1:]

    psf_OFF = open_fits(loadname % 'offaxis', verbose=False)
    print(f"Using off-axis PSF from {loadname%'offaxis'}")
    assert psf_OFF.ndim == 2, "off-axis PSF frame must be 2-dimensional"

    # crop everything to a common size
    min_crop = min(psf_OFF.shape[-1], disc_model.shape[-1], conf["ndet"])
    conf["ndet"] = min_crop
    if psf_OFF.shape[-1] > min_crop:
        psf_OFF = frame_crop(psf_OFF, min_crop, verbose=False)
    if disc_model.shape[-1] > min_crop:
        disc_model = frame_crop(disc_model, min_crop, verbose=False)

    # apply extinction if requested
    if extinction > 0 and source_xy is not None:
        source_x, source_y = source_xy
        extinction_factor = 10 ** (-0.4 * extinction)
        disc_model[source_y, source_x] *= extinction_factor
        print("Applied extinction of %.2f mag to source at (x, y) = (%d, %d)" % (extinction, source_x, source_y), flush=True)

    # load coronagraph transmission depending on the mode and band
    if conf["mode"] == "ELT" or conf["mode"] is None:  # no coronagraph
        transmission = None
    elif "VC" in conf["mode"]:
        if conf["band"] in ("L", "M"):
            conf['f_oat'] = conf['dir_input'] + 'optics/vc/oat_L_%s.fits' % conf['mode']
        elif conf["band"] == "N2":
            conf['f_oat'] = conf['dir_input'] + 'optics/vc/oat_N2_CVC.fits'
        transmission = open_fits(conf['f_oat'], verbose=False)
        print(f"Using transmission from {conf['f_oat']}")
    else: # TODO
        raise NotImplementedError(f"Mode {conf['mode']} not implemented for disc injection.")

    # open and crop sequence
    psf_ON = open_fits(loadname % 'onaxis', verbose=False)
    print(f"Using on-axis PSF from {loadname % 'onaxis'}")
    assert psf_ON.ndim == 3, "on-axis PSF cube must be 3-dimensional"

    if psf_ON.shape[-1] > min_crop:
        psf_ON = cube_crop_frames(psf_ON, min_crop, verbose=False)

    nframes = psf_ON.shape[0]

    # generate parallactic angles
    pa = paralang(npts=nframes, dec=conf["dec"], lat=conf["lat"], duration=int(nframes * conf["dit"]))

    # RDI handling
    if do_rdi:
        if rdi_duration is None:
            rdi_duration = int(nframes * conf["dit"] * 0.2)
            print("Note: rdi_duration was not set. Using 20% of the PSF sequence.", flush=True)

        n_ref_frames = int(round(rdi_duration / conf["dit"]))

        if rdi_mag == conf["mag"]:
            # always take the RDI block from the end of the sequence, leaving the science
            # sequence as one contiguous block from the start
            rdi_idx = slice(nframes - n_ref_frames, nframes)
            sci_idx = slice(0, nframes - n_ref_frames)

            psf_RDI = psf_ON[rdi_idx].copy()
            psf_ON = psf_ON[sci_idx]
            psf_OFF_rdi = psf_OFF.copy()
            pa = pa[sci_idx]

            print(f"RDI reference (mag={rdi_mag}) reused from the science cube: "
                  f"{n_ref_frames} frames removed from the end of the science sequence, "
                  f"no overlap with science frames.", flush=True)

        else:  # we load in a different PSF cube for the RDI reference and extract the requested duration
            rdi_grid_dir = os.path.join(
                db_path,
                f"{conf['band']}_{conf['mode']}_s=Q{seeing_q}_mag={rdi_mag}_{conf['duration']}s_{int(conf['dit'] * 1000)}ms"
            )
            rdi_loadname = os.path.join(rdi_grid_dir, '%s_PSF_bckg1_%s_%s.fits' % ('%s', conf['band'], conf['mode']))

            psf_OFF_rdi = open_fits(rdi_loadname % 'offaxis', verbose=False)
            print(f"Using RDI off-axis PSF from {rdi_loadname % 'offaxis'}")
            assert psf_OFF_rdi.ndim == 2, "RDI off-axis PSF frame must be 2-dimensional"

            psf_RDI = open_fits(rdi_loadname % 'onaxis', verbose=False)
            print(f"Using RDI on-axis PSF from {rdi_loadname % 'onaxis'}")
            assert psf_RDI.ndim == 3, "RDI on-axis PSF cube must be 3-dimensional"

            if psf_OFF_rdi.shape[-1] > min_crop:
                psf_OFF_rdi = frame_crop(psf_OFF_rdi, min_crop, verbose=False)
            if psf_RDI.shape[-1] > min_crop:
                psf_RDI = cube_crop_frames(psf_RDI, min_crop, verbose=False)

            start_idx = np.random.randint(low=0, high=psf_RDI.shape[0] - n_ref_frames + 1)
            psf_RDI = psf_RDI[start_idx:start_idx + n_ref_frames]

        if starphot is not None:
            _, _, ap_flux_rdi = psf_template(psf_OFF_rdi)
            psf_RDI *= starphot / ap_flux_rdi
        del psf_OFF_rdi

    # where the magic happens
    cube = cube_inject_fakedisk(
        disc_model,
        pa,
        transmission=transmission,
        psf=psf_OFF,
        normalize_psf=False,
        nproc=conf["cpu_count"],
        mask_val=0,
        edge_blend="interp",
        interp_zeros=True,
        ker=1,
        imlib=imlib,
    )

    # normalise by the stellar flux that we removed earlier (units of contrast)
    cube /= star_val

    # add the disc model
    psf_ON += cube
    del cube

    if starphot is not None:
        _, _, ap_flux = psf_template(psf_OFF)
        psf_ON *= starphot / ap_flux
    del psf_OFF

    if do_rdi:
        return psf_ON, pa, psf_RDI
    else:
        return psf_ON, pa


def _check_grid(mag : Union[float, int], band : str, mode : str) -> float:
    """
    Check if the provided magnitude is in the Liege grid. If not, round to the closest available magnitude. A float
    must be returned as the grid contains mag as a float in the filenames.

    Parameters
    ----------
    mag : float, int
        The magnitude to check.
    band : str
        The band of the observation (e.g., "L", "M", "N1", "N2".)
    mode : str
        The mode of the observation (e.g., "CVC", "RAVC".)

    Returns
    -------
    mag : float
        The closest available magnitude in the Liege grid for a given band and mode.
    """
    mag_og = float(mag)

    if band == "L" and mode == "CVC":
        allowed = np.arange(-1.5, 9.0 + 0.5, step=0.5)
    elif band == "L" and mode == "RAVC":
        allowed = np.arange(-1.5, 4.0 + 0.5, step=0.5)
    elif band == "M" and mode == "CVC":
        allowed = np.arange(-1.5, 7.0 + 0.5, step=0.5)
    elif band == "N1" and mode == "CVC":
        allowed = np.arange(-1.5, 4.0 + 0.5, step=0.5)
    elif band == "N2" and mode == "CVC":
        allowed = np.arange(-1.5, 3.0 + 0.5, step=0.5)
    else:
        raise ValueError(f"Band {band} and mode {mode} is not supported.")

    if mag not in allowed:
        mag = round(mag / 0.5) * 0.5

    # clamp to valid range
    if mag < allowed[0]:
        mag = allowed[0]
        print(f"Attention: mag={mag_og} is below minimum for {band}. Using mag={mag}")
    elif mag > allowed[-1]:
        mag  = allowed[-1]
        print(f"Attention: mag={mag_og} is above maximum for {band}. Using mag={mag}")
    elif mag != mag_og:
        print(f"Attention: mag={mag_og} rounded to nearest grid value mag={mag}")

    return float(mag)
