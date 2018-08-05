import collections
import os
import proper

conf = collections.OrderedDict()

conf['OUT_DIR'] = './output_files/'
conf['TMP_DIR'] = './temp_files/'
conf['INPUT_DIR'] = './input_files/'

conf['GRIDSIZE'] = 1024
conf['WAVELENGTH'] = 5.0*10**-6
conf['DIAM'] = 37.0
conf['R_OBSTR'] = 0.3
conf['SPIDERS_WIDTH'] = 0.60
conf['SPIDERS_ANGLE'] = [0, 60.0, 120.0]
conf['MIS_SEGMENTS_NU'] = 0
conf['PUPIL_FILE'] = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
conf['PREFIX'] = 'test'
conf['PIXEL_SCALE'] = 5.0
conf['F_LENS'] = 658.6
conf['DEBUG'] = False
conf['DEBUG_PRINT'] = False

conf['TILT_2D'] = [0.0, 0.0]
conf['TILT_CUBE'] = 'np.random.randn(10,2)'
conf['ATM_SCREEN_NO'] = 0.
conf['ATM_SCREEN_2D'] = 'metis_370P_35L_HCI_Feb18_rwf8160_cut.fits'
conf['ATM_SCREEN_CUBE'] = 'cube_atm_1000screens_Feb2018_RandomWind.fits'
conf['ISLAND_PISTON'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
conf['NCPA'] = ' 0.'

conf['CHARGE'] = 2

conf['PHASE_APODIZER_FILE'] = 'app_phase_cut.fits'
conf['APODIZER_MIS_ALIGN'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
conf['AMP_APODIZER'] = 0

conf['LYOT_STOP'] = True
conf['LS_PARA'] = [0.98, 0.03, 1.1]
conf['LS_MIS_ALIGN'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

conf['N_D'] = 512



collections.OrderedDict(sorted(conf.items()))


os.makedirs(conf['OUT_DIR'] , exist_ok=True)
os.makedirs(conf['TMP_DIR'], exist_ok=True)
os.makedirs(conf['INPUT_DIR'], exist_ok=True)


proper.print_it = False
