from astropy.io import fits 
import astropy.units as u
import time

# get header and data
case_name = 'drift' #all_effects_2048 #scao_only_2048
contrast = 'ADI'
mode = 'RAVC'
if contrast.upper() == 'ADI':
    contrast = contrast.upper()
elif contrast.lower() == 'raw':
    contrast = contrast.lower()
scao = 'compass3600s_samp300ms'    #'compass600s_samp100ms'  #'compass3600s_samp300ms'
band = 'M'
#band_bckg = '%s_bckg0'%band
band_bckg = '%s_bckg0'%band
#band_bckg = 'L_bckg1'
filename = 'cc_%s_dec-2.47_%s_%s_%s_%s.fits'%(scao, band_bckg, mode, case_name, contrast)

print(filename)
fits.open(filename)[0].header
hdu = fits.open(filename)[0]
hdr = hdu.header
data = hdu.data

# define new header elements
# hdr['bckg'] = False#True
# hdr['mag'] = None #if bckg is True else None
# hdr['date'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
# hdr['title'] = 'SCAO residuals only (1h, 300ms)'
# hdr['mode'] = 'CVC'
# hdr['band'] = 'N2a'#'N1a'#'N2'#'N1'#'M'#'L'
hdr['lam'] = 3.8128e-06#4.8091E-06#8.6016E-06#1.1236E-05#8.6469E-06#1.1187E-05  #1024 values
#hdr['lam'] = 1.1219E-05#3.8204e-06#4.7970E-06#8.6313E-06#1.1287E-05#8.6658E-06#  #2048 values
# hdr['diam'] = 37

# change header
if contrast == 'ADI':
    hdr['ylabel'] = '5-$\sigma$ sensitivity (contrast)'
elif contrast == 'raw':
    hdr['ylabel'] = 'Raw contrast'
hdr['xlabel'] ='Angular separation $[\lambda/D]$'
hdr2 = fits.Header({'xlabel':hdr['xlabel'], 'ylabel':hdr['ylabel'], \
        'title':hdr['title']})
hdr2['date'] = (hdr['date'], 'created')
hdr2['contrast'] = (contrast, 'contrast type: ADI, PCA, raw,...')
hdr2['mode'] = (mode, 'HCI mode')
hdr2['band'] = (hdr['band'], 'spectral band')
hdr2['lam'] = (hdr['lam'], 'wavelength in m')
hdr2['diam'] = (hdr['diam'], 'pupil diameter in m')
hdr2['bckg'] = (hdr['bckg'], 'shot noise (bool)')
hdr2['mag']  = (hdr['mag'], 'star magnitude (if bckg is True else None)')

#hdr2['title'] = 'static NCPA at HSF (>10c/p)'
#hdr2['title'] = 'quasi-static NCPA (0-0.01Hz) at HSF (>10c/p)'
#hdr2['title'] = 'quasi-static NCPA (0-0.01Hz) at LSF (<10c/p)'
#hdr2['title'] = 'quasi-static petal piston (0-0.01Hz) at HSF (<10c/p)'
#hdr2['title'] = 'quasi-static NCPA (<fc) at HSF (>10c/p)'
#hdr2['title'] = 'quasi-static NCPA (<fc) at LSF (<10c/p)'
#hdr2['title'] = 'dynamic NCPA (0.01-1Hz) at all SF'
#hdr2['title'] = 'dynamic petal piston (0.01-1Hz) at all SF'
#hdr2['title'] = 'pointing jitter (0.01-1Hz)'
#hdr2['title'] = 'pointing quasistatics (0-0.01Hz)'
#hdr2['title'] = 'SCAO-corrected static NCPA at LSF (<10c/p), 600s/100ms'
#hdr2['title'] = 'SCAO residuals only (1h, 300ms)'
#hdr2['title'] = '10min comparison w/ Faustine 20181008 (1h)'
#hdr2['title'] = 'various NCPA sources'
#hdr2['title'] = 'various error sources'
#hdr2['title'] = 'NCPA linear drift'
#hdr2['title'] = 'petal piston'
#hdr2['title'] = 'background limited perf. with SCAO residuals only'
#hdr2['title'] = 'background limited perf. with all error sources'
#hdr2['title'] = 'APP comparison'

#hdr2['title'] = 'CVC at %s band, SCAO residuals only'%band
#hdr2['title'] = 'CVC at %s band, with all effects'%band
#hdr2['title'] = 'CVC at %s band (geosnap), SCAO residuals only'%band
#hdr2['title'] = 'CVC at %s band (geosnap), with all effects'%band
#bandx = band

#hdr2['title'] = 'CVC at N2 band (aquarius), SCAO residuals only'
#hdr2['title'] = 'CVC at N2 band (aquarius), with all effects'
#bandx = 'N2'

hdr2['title'] = 'M-band RAVC - apodizer drift study'

hdr2['LABEL0'] = 'x-values'
hdr2['LABEL1'] = 'drift=1% ptv'
hdr2['LABEL2'] = 'drift=2% ptv'
hdr2['LABEL3'] = 'drift=3% ptv'

# hdr2['LABEL0'] = 'x-values'
# hdr2['LABEL1'] = '%s=3.5'%bandx
# hdr2['LABEL2'] = '%s=3.0'%bandx
# hdr2['LABEL3'] = '%s=2.5'%bandx
# hdr2['LABEL4'] = '%s=2.0'%bandx
# hdr2['LABEL5'] = '%s=1.5'%bandx
# hdr2['LABEL6'] = '%s=1.0'%bandx
# hdr2['LABEL7'] = '%s=0.5'%bandx
# hdr2['LABEL8'] = '%s=0.0'%bandx
# hdr2['LABEL9'] = '%s=-0.5'%bandx
# hdr2['LABEL10'] = '%s=-1.0'%bandx
# hdr2['LABEL11'] = '%s=-1.5'%bandx

#hdr2['LABEL1'] = 'ELT pupil'
#hdr2['LABEL2'] = 'binary pupil'

#hdr2['LABEL1'] = 'no background'
#hdr2['LABEL2'] = 'L = 5 (end-to-end)'
#hdr2['LABEL3'] = 'L = 7.5 (end-to-end)'
#hdr2['LABEL4'] = 'L = 10 (end-to-end)'
#hdr2['LABEL5'] = 'L = 5 (reconstructed)'
#hdr2['LABEL6'] = 'L = 7.5 (reconstructed)'
#hdr2['LABEL7'] = 'L = 10 (reconstructed)'


#hdr2['LABEL2'] = 'LSF (<10cpp) + LTF (<0.01Hz, 10nm rms)'
#hdr2['LABEL3'] = 'HSF (>10cpp) + LTF (<0.01Hz, 10nm rms)'
#hdr2['LABEL4'] = 'LSF (<10cpp) + HTF (>0.01Hz, 10nm rms)'
#hdr2['LABEL5'] = 'HSF (>10cpp) + HTF (>0.01Hz, 10nm rms)'


#hdr2['LABEL1'] = '0 nm rms'
#hdr2['LABEL2'] = '200 nm rms'
#hdr2['LABEL3'] = '125 nm rms'
#hdr2['LABEL4'] = '150 nm rms'
#hdr2['LABEL5'] = '175 nm rms'
#hdr2['LABEL6'] = '200 nm rms'
#hdr2['LABEL7'] = '225 nm rms'
#hdr2['LABEL8'] = '250 nm rms'

#hdr2['LABEL2'] = 'static NCPA at HSF, 35.9 nm rms'
#hdr2['LABEL3'] = 'quasi-static NCPA at LSF, 20 nm rms'
#hdr2['LABEL4'] = 'quasi-static NCPA at HSF, 20 nm rms'
#hdr2['LABEL5'] = 'dynamic NCPA at all SF, 40 nm rms'
#hdr2['LABEL6'] = 'all NCPA effects combined'

# hdr2['LABEL2'] = 'pointing jitter 2 mas rms'
# hdr2['LABEL3'] = 'pointing quasistatics 0.4 mas rms'
# hdr2['LABEL4'] = 'pupil drift (RAVC) 2% ptv'
# hdr2['LABEL5'] = 'NCPA + petal piston (see box)'
# hdr2['LABEL6'] = 'ALL EFFECTS'
# hdr2['LABEL7'] = 'ALL EFFECTS + 7 misaligned segments'

## edit data
data2 = data[0:1,:]
#data2 = np.vstack((data2, data[-2,:]))

#cube0 = fits.getdata('cc_compass3600s_samp300ms_dec-2.47_L6_bckg0_RAVC_SCAOonly_raw.fits')
#cube1 = fits.getdata('cc_raw_L_RAVC_misaligned_segments_7_flower.fits')
#cube2 = fits.getdata('cc_raw_L_RAVC_pointing_jitter_0.01-1Hz_2.fits')
#cube3 = fits.getdata('cc_raw_L_RAVC_pointing_quasistatics_0-0.01Hz_0.4.fits')
#cube4 = fits.getdata('cc_raw_L_RAVC_pupil_drift_RAVC_2.fits')
##cube5 = fits.getdata('cc_raw_L_RAVC_ncpa+petal_piston_20.fits')
#cube5 = fits.getdata('cc_raw_L_RAVC_petal_piston_drift_10.fits')
##cube6 = fits.getdata('cc_raw_L_RAVC_ncpa_lin_drift_ptv_20.fits')
#cube6 = fits.getdata('cc_raw_L_RAVC_all_ncpa_20191011.fits')
#cube7 = fits.getdata('cc_raw_L_RAVC_ALL_EFFECTS_20191017.fits')

#cube0 = fits.getdata('cc_compass3600s_samp300ms_dec-2.47_L6_bckg0_RAVC_SCAOonly_ADI.fits')
#cube1 = fits.getdata('cc_compass1h_samp300ms_dec-2.47_L_mag6_bckg0_RAVC_all_effects_2.fits')
#cube2 = fits.getdata('cc_compass3600s_samp300ms_ADI3600s_samp300ms_avg0ms_dec-2.47deg_L_mag6_bckg0_RAVC_71.fits')
#cube3 = fits.getdata('cc_compass3600s_samp300ms_dec-2.47_L6_bckg0_RAVC_pointing_quasistatics_0-0.01Hz_ADI.fits')
#cube4 = fits.getdata('cc_compass1h_samp300ms_dec-2.47_L_mag6_bckg0_RAVC_all_effects_5.fits')
##cube5 = fits.getdata('cc_raw_L_RAVC_ncpa+petal_piston_20.fits')
#cube5 = fits.getdata('cc_compass1h_samp300ms_dec-2.47_L_mag6_bckg0_RAVC_all_effects_6.fits')
##cube6 = fits.getdata('cc_raw_L_RAVC_ncpa_lin_drift_ptv_20.fits')
#cube6 = fits.getdata('cc_compass3600s_samp300ms_dec-2.47_L6_bckg0_RAVC_20191011_ADI.fits')
#cube7 = fits.getdata('cc_compass3600s_samp300ms_ADI3600s_samp300ms_avg0ms_dec-2.47deg_L_mag6_bckg0_RAVC_78.fits')

if contrast == 'ADI':
    update_x = True
#    for i in np.arange(3.5,-2,-.5):
#        cube = fits.getdata('cc_%s_dec-2.47_%smag%s_bckg1_CVC_%s.fits'%(scao,band,round(i,1),case_name)).T
    for i in np.arange(1,4,1):
        cube = fits.getdata('cc_%s_dec-2.47_%smag6_bckg0_RAVC_%s_%s.fits'%(scao,band,case_name,i)).T
        # X-axis in lam/D
        if update_x is True:
            lamD2asec = hdr['lam']/hdr['diam']*u.rad.to('arcsec') # lam/D to arcsec
            data2 = cube[4,:]/lamD2asec
            update_x = False
        data2 = np.vstack((data2, cube[1,:]))
    # cube0 = fits.getdata('cc_%s_ADI3600s_samp300ms_avg0ms_dec-2.47deg_%s_mag3.5_bckg1_CVC_cvc_%s.fits'%(scao,band,case_name)).T
    # cube1 = fits.getdata('cc_%s_ADI3600s_samp300ms_avg0ms_dec-2.47deg_%s_mag3_bckg1_CVC_cvc_%s.fits'%(scao,band,case_name)).T
    # cube2 = fits.getdata('cc_%s_ADI3600s_samp300ms_avg0ms_dec-2.47deg_%s_mag2.5_bckg1_CVC_cvc_%s.fits'%(scao,band,case_name)).T
    pass
elif contrast == 'raw':
    cube = fits.getdata('cc_raw_L_RAVC_SCAO_ONLY_0.fits')
    cube0 = fits.getdata('cc_raw_L_RAVC_ncpa+petal_piston_20191205.fits')
    cube1 = fits.getdata('cc_raw_L_RAVC_ALL_EFFECTS_20191205.fits')
    cube2 = fits.getdata('cc_raw_L_RAVC_ALL_EFFECTS+MIS_SEG_20191205.fits')
    cube3 = fits.getdata('cc_raw_L_RAVC_pointing_quasistatics_0-0.01Hz_0.4.fits')
    cube4 = fits.getdata('cc_raw_L_RAVC_pointing_jitter_0.01-1Hz_2.fits')
    cube5 = fits.getdata('cc_raw_L_RAVC_pupil_drift_RAVC_2.fits')
    pass



if False:
    data2 = np.vstack((data2, cube0[1,:]))
    data2 = np.vstack((data2, cube1[1,:]))
    data2 = np.vstack((data2, cube2[1,:]))
    data2 = np.vstack((data2, cube3[1,:]))
    data2 = np.vstack((data2, cube4[1,:]))
    data2 = np.vstack((data2, cube5[1,:]))
    data2 = np.vstack((data2, cube6[1,:]))
    data2 = np.vstack((data2, cube7[1,:]))
    data2 = np.vstack((data2, cube8[1,:]))

#data2 = np.vstack((data2, cube1[1,:]))
#data2 = np.vstack((data2, cube1[1,:]))
#data2 = np.vstack((data2, cube2[1,:]))

#C = 9.4e-09 # [median(mag=10) - median(no bckg)] / 10^(0.4*mag)
#data2 = np.vstack((data2, data2[1,:] + C*10**(0.4*5)))
#data2 = np.vstack((data2, data2[1,:] + C*10**(0.4*7.5)))
#data2 = np.vstack((data2, data2[1,:] + C*10**(0.4*10)))


#data2 = np.vstack((data2, data2[1,:]*C*10**(0.4*7.5)))
#data2 = np.vstack((data2, data2[1,:]*C*10**(0.4*10)))


# save fits file
hdu = fits.PrimaryHDU(data2, header=hdr2) # select data2 / hdr2
hdu.writeto(filename, overwrite=True)
print('new header')
fits.open(filename)[0].header

