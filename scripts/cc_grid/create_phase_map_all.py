from heeps.util.multiCPU import multiCPU
from heeps.util.freq_decomp import fit_zer, remove_zernike
from heeps.util.img_processing import resize_cube
import numpy as np
import os
from astropy.io import fits
from copy import deepcopy
import proper

# Temporal sampling of the SCAO phase cubes [s]
# Used to convert NCPA loop frequency [Hz] → frame-decimation factor
PHASE_SCREEN_DIT = 0.1   # s  (i.e. 10 Hz cubes)


class PhaseCubeGenerator():
    """
    Build and write a combined NCPA+SCAO+WV phase cube to a FITS file.

    All selection logic (filenames, noise levels, loop parameters) is resolved
    upstream by Sacred and passed in via the config dict.  This class is
    responsible only for the numerical assembly and caching of intermediate
    products (Zernike polynomial fits).

    Typical usage (from Sacred run_simulation)::

        sigLF     = get_wfe_tip_tilt(band)
        sigHF     = get_wfe_higher_order(band, magnitude, mode, ncpa_freq)
        generator = PhaseCubeGenerator(_config)
        generator.run(sigLF, sigHF)
    """

    # RMS of the reference SCAO phase screen used for WV normalisation [nm]
    # (see temporal_rms_matisse.ipynb)
    TEMPORAL_RMS_REF = 8814.1

    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict
            Sacred ``_config`` dictionary.  Required keys (all provided by
            ``derived_config`` in the Sacred experiment):
            ``f_scao``, ``f_wv``, ``f_cbw``, ``f_pupil``, ``f_phase``,
            ``band``, ``npupil``, ``wv_rms``, ``ncpa_freq``,
            ``ncpa`` (sub-dict with keys ``lag``, ``nmodes``, ``gain_I``).
        """
        # Input / output filenames
        dir_inputs = config['dir_input'] + '/'
        self.scao_filename  = dir_inputs + config['f_scao']
        self.wv_filename    = dir_inputs + config['f_wv']
        self.cbw_filename   = dir_inputs + config['f_cbw']
        self.pup_filename   = dir_inputs + config['f_pupil']
        self.phase_filename = dir_inputs + config['f_phase']

        # Scalar parameters from config
        self.band     = config['band']
        self.npupil   = config['npupil']
        self.wv_rms   = config['wv_rms']

        # NCPA closed-loop parameters (from ncpa ingredient)
        # nzern = number of ALF-corrected higher-order modes (excludes piston, tip, tilt)
        # total Zernike array size = 3 + nzern  (index 0=piston, 1-2=tip-tilt, 3…=HO)
        self.nzern    = config['ncpa']['nmodes']
        self.lag      = config['ncpa']['lag']
        self.gain_I   = config['ncpa']['gain_I']

        # Frame-decimation factors: how many phase-screen frames per correction step.
        # Phase screens are sampled at PHASE_SCREEN_DIT; ncpa_freq is the loop rate.
        # Both tip-tilt (LF) and higher-order (HF) use the same NCPA frequency here;
        # extend with a separate key if they need to differ.
        self.ncpa_freq = config['ncpa_freq']
        self.n_decim   = max(1, round(1.0 / (PHASE_SCREEN_DIT * self.ncpa_freq)))

        self.cpu_count=config['cpu_count']

    # ------------------------------------------------------------------
    def run(self, sigLF, sigHF):
        """
        Build the combined phase cube and write it to ``self.phase_filename``.

        If the output file already exists it is left untouched (cache hit).

        Parameters
        ----------
        sigLF : float
            Tip-tilt sensor noise [nm rms]  (from ``get_wfe_tip_tilt``).
        sigHF : float
            Higher-order sensor noise [nm rms]  (from ``get_wfe_higher_order``).

        Returns
        -------
        str
            Path of the (new or existing) output FITS file.
        """
        # sigLF and sigHF need to be in units of meter rms for load_zpols_integ
        sigLF *= 1e-9
        sigHF *= 1e-9
        if os.path.isfile(self.phase_filename):
            print('file already exists: ' + self.phase_filename)
            return self.phase_filename

        print('write to ' + self.phase_filename)
        pup, ncpa = self.build_all_phase()

        zpols_integ = self.load_zpols_integ(pup, ncpa, sigLF, sigHF)

        wf = proper.prop_begin(1, 1, self.npupil, 1)  # initial wavefront

        print(f'[DEBUG] {pup.shape}, {ncpa.shape}, {len(zpols_integ)}')

        _, HSF = multiCPU(remove_zernike, multi_out=True, verbose=True,
            posargs=[deepcopy(wf), pup],
            posvars=[ncpa, zpols_integ],
            cpu_count=self.cpu_count,
            reuse_pool=False)

        hdr = self._build_fits_header(sigLF, sigHF)
        fits.writeto(self.phase_filename, np.float32(HSF), hdr)

        return self.phase_filename


    # ------------------------------------------------------------------
    def build_all_phase(self):
        """
        Load and combine SCAO, water-vapour, and CBW-NCPA phase cubes.

        The WV cube is scaled so that its temporal RMS matches ``self.wv_rms``
        (the value from the Sacred config).

        Returns
        -------
        pup : ndarray
            Binary pupil mask (values 0 or 1).
        ncpa : ndarray
            Combined phase cube [same units as input FITS files].
        """
        scao       = fits.getdata(self.scao_filename)
        wv_unscaled = fits.getdata(self.wv_filename)
        cbw        = fits.getdata(self.cbw_filename)

        if cbw.shape[0] > scao.shape[0]:
            nframes = cbw.shape[0] - scao.shape[0]
            print(f'Warning: trimming cbw ncpa by {nframes} frames')
            cbw = cbw[:-nframes]

        pup = fits.getdata(self.pup_filename)
        pup[pup <= 0.5] = 0

        wv   = wv_unscaled * self.wv_rms / self.TEMPORAL_RMS_REF
        ncpa = scao + wv + cbw

        return pup, ncpa


    # ------------------------------------------------------------------
    def _build_fits_header(self, sigLF, sigHF):
        """
        Build a FITS header recording the key simulation parameters.

        Parameters
        ----------
        sigLF : float
            Tip-tilt sensor noise (converted to [m rms] in ``run``).
        sigHF : float
            Higher-order sensor noise (converted to [m rms] in ``run``).
        """
        hdr = fits.Header()
        entries = {
            'BAND':    (self.band,                          'observing band'),
            'NPUPIL':  (self.npupil,                        'pupil grid size [pix]'),
            'WVRMS':   (self.wv_rms,                        'water-vapour OPD rms [nm]'),
            'WVRMSREF':(self.TEMPORAL_RMS_REF,              'ref. temporal rms for WV scaling [nm]'),
            'NZERN':   (self.nzern,                         'ALF HO Zernike modes (no piston/tip/tilt)'),
            'LAG':     (self.lag,                           'NCPA loop lag [frames]'),
            'GAIN_I':  (self.gain_I,                        'integrator gain'),
            'NCPAFREQ':(self.ncpa_freq,                     'NCPA loop frequency [Hz]'),
            'NDECIM':  (self.n_decim,                       'decim. factor (phase frames per corr step)'),
            'SIGLF':   (np.round(sigLF*1e9, decimals=1),        'tip-tilt sensor noise [nm rms]'),
            'SIGHF':   (np.round(sigHF*1e9, decimals=1),        'higher-order sensor noise [nm rms]'),
            'F_SCAO':  (os.path.basename(self.scao_filename), 'SCAO phase cube filename'),
            'F_WV':    (os.path.basename(self.wv_filename),   'water-vapour phase cube filename'),
            'F_CBW':   (os.path.basename(self.cbw_filename),  'CBW NCPA phase cube filename'),
        }
        for key, (value, comment) in entries.items():
            hdr.set(key, value, comment)
        return hdr


    # ------------------------------------------------------------------
    def load_zpols(self, pup, ncpa):
        """
        Fit Zernike modes to every frame of the NCPA cube.

        Results are cached to a sidecar ``.fits`` file next to
        ``self.phase_filename`` so that repeated runs are fast.

        Parameters
        ----------
        pup : ndarray
            Binary pupil mask.
        ncpa : ndarray
            Combined phase cube.

        Returns
        -------
        zpols : ndarray, shape (nframes, nzern)
        """
        zpols_name = (self.phase_filename.replace('.fits', '')
                      + f'_zpols_{self.nzern}.fits')

        if os.path.isfile(zpols_name):
            print('    getdata ' + zpols_name)
            return fits.getdata(zpols_name)

        print('    writeto ' + zpols_name)
        zpols = multiCPU(fit_zer,
                         posargs=[pup, self.npupil / 2, 3 + self.nzern],
                         posvars=[ncpa],
                         case='get zpols',
                         cpu_count=self.cpu_count,
                         reuse_pool=False)
        fits.writeto(zpols_name, np.float32(zpols))
        return zpols


    # ------------------------------------------------------------------
    def load_zpols_integ(self, pup, ncpa, sigLF, sigHF):
        """
        Emulate a focal-plane closed-loop NCPA correction with an integrator.

        The correction is applied separately for tip-tilt (modes 1-2) and
        higher-order modes (3 … nzern-1) using the same frame-decimation factor
        ``self.n_decim`` and lag ``self.lag``.  Sensor noise is injected as
        zero-mean Gaussian noise with standard deviations ``sigLF`` and ``sigHF``.

        Results are cached to a sidecar ``.fits`` file.

        Parameters
        ----------
        pup   : ndarray  – binary pupil mask
        ncpa  : ndarray  – combined phase cube
        sigLF : float    – tip-tilt sensor noise [m rms]
        sigHF : float    – higher-order sensor noise [m rms]

        Returns
        -------
        zpols_integ : ndarray, shape (nframes, nzern)
        """
        nLF = nHF = self.n_decim
        zpols_integ_name = (
            self.phase_filename.replace('.fits', '')
            + f'_zpols_{self.nzern}_nLF_{nLF}_nHF_{nHF}'
            + f'_G_{self.gain_I}_lag_{self.lag}'
            + f'_sigLF_{sigLF:.1f}_sigHF_{sigHF:.1f}.fits'
        )

        if os.path.isfile(zpols_integ_name):
            print('  getdata ' + zpols_integ_name)
            return fits.getdata(zpols_integ_name)

        print('  write to ' + zpols_integ_name)
        zpols       = self.load_zpols(pup, ncpa)
        zpols_integ = np.zeros(zpols.shape)
        nframes     = len(zpols)

        # Piston: pass through uncorrected
        zpols_integ[:, 0] = zpols[:, 0]

        # Closed-loop correction: tip-tilt (modes 1-2) and higher order (3…3+nzern-1)
        for m, freq, sig in zip(
                [range(1, 3), range(3, 3 + self.nzern)],
                [nLF,         nHF],
                [sigLF,       sigHF]):
            for n in range(freq + self.lag, nframes, freq):
                error = (
                    np.mean(zpols[n-freq-self.lag : n-self.lag, m]
                            - zpols_integ[n-freq-self.lag : n-self.lag, m], axis=0)
                    + np.random.normal(0, sig, (1, len(m)))
                )
                zpols_integ[n : n+freq, m] = zpols_integ[n-1, m] + self.gain_I * error

        fits.writeto(zpols_integ_name, np.float32(zpols_integ))
        return zpols_integ

