from heeps.util.img_processing import crop_img, resize_img
from heeps.util.multiCPU import multiCPU
from heeps.util.notify import notify
import os
import numpy as np
from astropy.io import fits

# email when finished
conf = dict(send_to = 'cdelacroix@uliege.be',
            send_message = 'cube calculation finished OK.')

# useful inputs
os.chdir(os.path.normpath(os.path.expandvars('$HOME/heeps_metis/input_files')))
output_folder = 'wavefront/dfull'
prefix = 'Residual_phase_screen_'#'tarPhase_1hr_100ms_'
suffix = 'ms'
duration = 600#3600#
samp = 100#300#
start = 4101#101#2101#
nimg = 720
band = 'L'
npupil = 285#720#
pad_frame = False
cpu_count = None

input_folders = [
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag7.0_2024-03-11T20:52:15/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag7.0_2024-03-11T20:57:01/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag7.5_2024-03-11T22:26:43/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag7.5_2024-03-11T22:27:14/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag8.0_2024-03-12T09:07:29/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag8.0_2024-03-12T11:19:53/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag8.5_2024-03-12T11:45:58/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag8.5_2024-03-26T16:56:08/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag7.0_2024-03-12T12:03:32/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag7.0_2024-03-12T14:23:37/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag7.5_2024-03-12T14:24:27/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag7.5_2024-03-12T15:33:13/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag8.0_2024-03-12T20:17:40/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag8.0_2024-03-13T00:23:49/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag8.5_2024-03-13T02:20:02/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ2_z30d_Kmag8.5_2024-03-13T04:25:53/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag7.0_2024-03-13T04:46:30/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag7.0_2024-03-13T07:52:16/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag7.5_2024-03-13T07:55:49/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag7.5_2024-03-13T08:33:05/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag8.0_2024-03-13T13:18:56/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag8.0_2024-03-13T14:16:52/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag8.5_2024-03-13T19:26:30/residualPhaseScreens',
    '/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_sQ3_z30d_Kmag8.5_2024-03-13T21:17:10/residualPhaseScreens',
]

tags = [
    'Dfull_20240313_default_D_smedian_z30d_Kmag7.0_v1',
    'Dfull_20240313_default_D_smedian_z30d_Kmag7.0_v2',
    'Dfull_20240313_default_D_smedian_z30d_Kmag7.5_v1',
    'Dfull_20240313_default_D_smedian_z30d_Kmag7.5_v2',
    'Dfull_20240313_default_D_smedian_z30d_Kmag8.0_v1',
    'Dfull_20240313_default_D_smedian_z30d_Kmag8.0_v2',
    'Dfull_20240313_default_D_smedian_z30d_Kmag8.5_v1',
    'Dfull_20240313_default_D_smedian_z30d_Kmag8.5_v2',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag7.0_v1',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag7.0_v2',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag7.5_v1',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag7.5_v2',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag8.0_v1',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag8.0_v2',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag8.5_v1',
    'Dfull_20240313_default_D_sQ2_z30d_Kmag8.5_v2',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag7.0_v1',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag7.0_v2',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag7.5_v1',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag7.5_v2',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag8.0_v1',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag8.0_v2',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag8.5_v1',
    'Dfull_20240313_default_D_sQ3_z30d_Kmag8.5_v2',
]

#input_folders = ['/mnt/disk12tb/METIS/METIS_COMPASS_RAW_PRODUCTS/Hinkley_032024/gorban_metis_Dfull_default_D_smedian_z30d_Kmag8.5_2024-03-26T16:56:08/residualPhaseScreens']
#tags = ['Dfull_20240313_default_D_smedian_z30d_Kmag8.5_v2']


for input_folder,tag in zip(input_folders,tags):

    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cbasic_2020-10-16T10:25:14/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cbasic_2020-11-05T12:40:27/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cbasic_2020-11-30T20:52:24/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cbasic_uncorrected_2021-06-01T12:02:36/residualPhaseScreens'
    #input_folder = '/mnt/disk12tb/Users/gorban/METIS/METIS_COMPASS/gorban_metis_baseline_Cfull_noWtt_2021-10-07T09:00:32/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cfull_noWtt_2021-10-04T12:28:42/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cfull_K10mag_2022-01-24T11:47:12/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_baseline_Cfull_phasescreen_only_p1..6_2022-02-25T14:03:22/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_Cfull_1hr_2022-05-12T02:19:28/residualPhaseScreens'
    #input_folder = '/mnt/disk4tb/METIS/METIS_COMPASS_RAW_PRODUCTS/gorban_metis_Cfull_FDR_1hr_2022-09-29T12:46:22/residualPhaseScreens'

    #tag = 'Dfull_20220929'#'Cfull_20220929'#'Cfull_20220512'#'Cfull_20220225'#'Cfull_20211004'#'Cfull_20211007'#'Cbasic_20201130'
    savename = 'cube_%s_%ss_%sms_scao_only_%s_%s.fits'%(tag, duration, samp, band, npupil)
    #savename = 'cube_%s_%ss_%sms_0piston_meters_scao_only_%s_WV.fits'%(tag, duration, samp, npupil)

    # mask
    mask = fits.getdata(os.path.join(input_folder, 'Telescope_Pupil.fits'))
    mask = crop_img(mask, nimg)
    mask_pupil = np.rot90(resize_img(mask, npupil))
    fits.writeto(os.path.join(output_folder, 'mask_%s_%s_%s.fits'%(tag, band, npupil)), np.float32(mask_pupil), overwrite=True)

    # filenames
    nframes = len([name for name in os.listdir(input_folder) if name.startswith(prefix)])
    #nframes = 12000 # 6000
    print('%s frames'%nframes)
    frames = [str(frame).zfill(6) if pad_frame is True else str(frame) \
        for frame in range(start, start + nframes*samp, samp)]
    filenames = np.array([os.path.join(input_folder, '%s%s%s.fits'%(prefix, frame, suffix)) \
        for frame in frames])

    def remove_piston(filename):
        data = np.float32(fits.getdata(filename))
        data = crop_img(data, nimg)
        data -= np.mean(data[mask!=0])      # remove piston
        data[mask==0] = 0
        data = resize_img(data, npupil)
        data = np.rot90(data) * 1e-6        # rotate, convert to meters
        return data

    # create cube
    cube = multiCPU(remove_piston, posvars=[filenames], case='create cube SCAO', cpu_count=cpu_count)
    print(cube.shape)
    fits.writeto(os.path.join(output_folder, savename), np.float32(cube), overwrite=True)
    notify(conf['send_message'], conf['send_to'])