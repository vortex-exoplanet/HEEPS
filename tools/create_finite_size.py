import heeps
from heeps.util.img_processing import get_radial_profile
from heeps.util.coord import disk_coord, mas2rms
import numpy as np
from astropy.io import fits


def raw_contrast(psf_OFF, psf_ON, ndet=403, pscale=5.47, xbin=1):
    rim = int((ndet - 1)/2) # image radius
    (xo,yo) = (rim, rim)
    xval = np.arange(rim)*pscale
    # on-axis PSFs (cube)
    if psf_ON.ndim != 3:
        psf_ON = np.array(psf_ON, ndmin=3)
    # average
    psf_ON_avg = np.mean(psf_ON, 0)
    # radial profiles
    yoff = get_radial_profile(psf_OFF, (xo,yo), xbin)[:-1]
    yon = get_radial_profile(psf_ON_avg, (xo,yo), xbin)[:-1]
    # normalize by the peak of the off-axis PSF
    yval = yon / np.max(yoff)
    return np.array([xval, yval]), psf_ON_avg


''' IMG mode, no error '''
print('''\nIMG mode, no error ''')

conf = dict(
    cpu_count = None,
    #dir_current = '$HOME/INSTRUMENTS/METIS/heeps_analysis',
    #dir_current = '/mnt/disk4tb/METIS/heeps_analysis',
    dir_output = 'output_files/create_finite_size',
    file_pupil = 'pupils/ELT_fullM1.fits',
    file_lyot_stop = 'pupils/ls_ravc_circ_285.fits',
    add_phase = False,
    band = 'L',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
psf = heeps.optics.detector(wf, savefits=False, verbose=True, **conf)
IMG_no_error, psf_ON_avg = raw_contrast(psf, psf)
fits.writeto('psf_IMG_no_error.fits', np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_IMG_no_error.fits', IMG_no_error, overwrite= True)

''' HCI point-like '''
print('''\nHCI point-like ''')

# no error

mode = 'CLC'
lyotstop = 'ls_clc_allglass_285.fits'

conf = dict(
    cpu_count = None,
    #dir_current = '$HOME/INSTRUMENTS/METIS/heeps_analysis',
    #dir_current = '/mnt/disk4tb/METIS/heeps_analysis',
    dir_output = 'output_files/create_finite_size',
    file_pupil = 'pupils/ELT_fullM1.fits',
    file_lyot_stop = 'pupils/%s'%lyotstop,
    add_phase = False,
    band = 'L',
    mode = mode,
    clc_diam = 4,
    ravc_calc = False,
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psf = heeps.wavefront.propagate_one(wf, onaxis=True, savefits=True, verbose=True, **conf)
HCI_PS_no_error, psf_ON_avg = raw_contrast(off, psf)
fits.writeto('psf_%s_PS_no_error.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PS_no_error.fits'%mode, HCI_PS_no_error, overwrite= True)

# SCAO only
print('\nscao only')

conf.update(
    nframes = 100,
    add_phase = True,
    file_phase = 'WFerrors/cube_Cbasic_20201130_3600s_300ms_0piston_meters_scao_only_285.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_PS_scao, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_PS_scao.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PS_scao.fits'%mode, HCI_PS_scao, overwrite= True)

# SCAO + pointing
print('\npointing')

conf.update(
    add_point_err = True,
    file_point_err = 'WFerrors/point_all_3600s_300ms.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_PS_point, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_PS_point.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PS_point.fits'%mode, HCI_PS_point, overwrite= True)

''' HCI, alpha Cen, r = 4.26 mas '''
print('''\nHCI, alpha Cen, r = 4.26 mas ''')

# no error

conf = dict(
    cpu_count = None,
    #dir_current = '$HOME/INSTRUMENTS/METIS/heeps_analysis',
    #dir_current = '/mnt/disk4tb/METIS/heeps_analysis',
    dir_output = 'output_files/create_finite_size',
    file_pupil = 'pupils/ELT_fullM1.fits',
    file_lyot_stop = 'pupils/%s'%lyotstop,
    add_phase = False,
    band = 'L',
    mode = mode,
    clc_diam = 4,
    ravc_calc = False,
    fp_offsets = mas2rms(disk_coord(4.26, nr=4), 39.9988)
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psf = heeps.wavefront.propagate_one(wf, onaxis=True, savefits=True, verbose=True, **conf)
HCI_AC_no_error, psf_ON_avg = raw_contrast(off, psf)
fits.writeto('psf_%s_AC_no_error.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_AC_no_error.fits'%mode, HCI_AC_no_error, overwrite= True)

# SCAO only
print('\nscao only')

conf.update(
    nframes = 100,
    add_phase = True,
    file_phase = 'WFerrors/cube_Cbasic_20201130_3600s_300ms_0piston_meters_scao_only_285.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_AC_scao, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_AC_scao.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_AC_scao.fits'%mode, HCI_AC_scao, overwrite= True)

# SCAO + pointing
print('\npointing')

conf.update(
    add_point_err = True,
    file_point_err = 'WFerrors/point_all_3600s_300ms.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_AC_point, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_AC_point.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_AC_point.fits'%mode, HCI_AC_point, overwrite= True)

''' HCI, pi1 Gru, r = 9.18 mas '''
print('''\nHCI, pi1 Gru, r = 9.18 mas ''')

# no error

conf = dict(
    cpu_count = None,
    #dir_current = '$HOME/INSTRUMENTS/METIS/heeps_analysis',
    #dir_current = '/mnt/disk4tb/METIS/heeps_analysis',
    dir_output = 'output_files/create_finite_size',
    file_pupil = 'pupils/ELT_fullM1.fits',
    file_lyot_stop = 'pupils/%s'%lyotstop,
    add_phase = False,
    band = 'L',
    mode = mode,
    clc_diam = 4,
    ravc_calc = False,
    fp_offsets = mas2rms(disk_coord(9.18, nr=4), 39.9988)
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psf = heeps.wavefront.propagate_one(wf, onaxis=True, savefits=True, verbose=True, **conf)
HCI_PG_no_error, psf_ON_avg = raw_contrast(off, psf)
fits.writeto('psf_%s_PG_no_error.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PG_no_error.fits'%mode, HCI_PG_no_error, overwrite= True)

# SCAO only
print('\nscao only')

conf.update(
    nframes = 100,
    add_phase = True,
    file_phase = 'WFerrors/cube_Cbasic_20201130_3600s_300ms_0piston_meters_scao_only_285.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_PG_scao, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_PG_scao.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PG_scao.fits'%mode, HCI_PG_scao, overwrite= True)

# SCAO + pointing
print('\npointing')

conf.update(
    add_point_err = True,
    file_point_err = 'WFerrors/point_all_3600s_300ms.fits',
)
conf = heeps.config.read_config(verbose=False, **conf)
conf = heeps.config.update_config(saveconf=False, verbose=True, **conf)
phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)
wf = heeps.pupil.pupil(savefits=False, verbose=True, **conf)
off = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)
psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
    amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
    savefits=True, verbose=True, **conf)
HCI_PG_point, psf_ON_avg = raw_contrast(off, psfs)
fits.writeto('psf_%s_PG_point.fits'%mode, np.float32(psf_ON_avg), overwrite= True)
fits.writeto('cc_raw_%s_PG_point.fits'%mode, HCI_PG_point, overwrite= True)



