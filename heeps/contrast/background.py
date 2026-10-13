import numpy as np
from astropy.io import fits
from heeps.contrast.sim_heeps import sim_heeps

def background(psf_ON, psf_OFF, header=None, mode='RAVC', lam=3.8e-6, dit=0.3, 
        mag=5, mag_ref=0, flux_star=9e10, flux_bckg=9e4, app_single_psf=0.48, 
        f_vc_trans=None, f_app_trans=None, seed=123456, verbose=False, 
        call_ScopeSim=False, **conf):

    """
    This function applies background and photon noise to intup PSFs (off-axis
    and on-axis), incuding transmission, star flux, and components transmittance.

    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        mode (str):
            HCI mode: RAVC, CVC, APP, CLC
        lam (float):
            wavelength in m
        dit (float):
            detector integration time in s
        mag (float):
            star magnitude
        mag_ref (float):
            reference magnitude for star and background fluxes
        flux_star (float):
            star flux at reference magnitude
        flux_bckg (float):
            background flux at reference magnitude
        app_single_psf (float):
            APP single PSF (4% leakage)
        f_vc_trans (str):
            path to VC transmittance fits file
        f_app_trans
            path to APP transmittance fits file
        seed (int):
            seed used by numpy.random process
        call_ScopeSim (bool):
            true if interfacing ScopeSim

    Return:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
    """

    # calculate offaxis-PSF transmission
    thruput = np.sum(psf_OFF)
    # load mask transmittance
    if 'VC' in mode:
        data = fits.getdata(f_vc_trans)
        mask_trans = np.interp(lam*1e6, data[0], data[1])
    elif 'APP' in mode:
        data = fits.getdata(f_app_trans)
        mask_trans = np.interp(lam*1e6, data[0], data[1])
    else:
        mask_trans = 1
    # apply correction for APP single PSF (~48%)
    if 'APP' in mode:
        psf_OFF *= app_single_psf
        psf_ON *= app_single_psf

    # scopesim-heeps interface
    if call_ScopeSim is True:

        hdu=fits.PrimaryHDU(data=psf_ON)
        hdu.header.append("CDELT1")
        hdu.header.append("CDELT2")
        hdu.header.append("CTYPE1")
        hdu.header.append("CTYPE2")
        hdu.header.append("CUNIT1")
        hdu.header.append("CUNIT2")
        hdu.header.append("CPIX1")
        hdu.header.append("CPIX2")
        hdu.header.append("CVAL1")
        hdu.header.append("CVAL2")
        hdu.header.append("BUNIT")
        hdu.header.append("OFFAXISTRANS")
        hdu.header.append("MODE")
        hdu.header.append("BAND")
        hdu.header.append("LAM")
        hdu.header.append("FILTER")
        if conf['ScopeSim_LMS'] == True:
            hdu.header.append("NAXIS3")
            hdu.header.append("CDELT3")
            hdu.header.append("CUNIT3")
            hdu.header.append("CPIX3")
            hdu.header.append("CVAL3")

        for band in conf['band_specs']:
            print(band,conf['band_specs'][band]['lam'])
            if np.abs(conf['band_specs'][band]['lam'] - lam)<0.05e-6:
                hdu.header['BAND']=band
                hdu.header['LAM']=lam
                
                
        if hdu.header['BAND']=="L":
            hdu.header["FILTER"]="HCI_L_long"
        elif hdu.header['BAND']=="M":
            hdu.header["FILTER"]="CO_ref"
        elif hdu.header['BAND']=="N1":
            hdu.header["FILTER"]="N1"
        elif hdu.header['BAND']=="N2":
            hdu.header["FILTER"]="N2"
        elif hdu.header['BAND']=="N1a":
            hdu.header["FILTER"]="N1" # no N2a exists in Scopesim
        elif hdu.header['BAND']=="N2a":
            hdu.header["FILTER"]="N2" # no N2a exists in Scopesim
        hdu.header['MODE']=mode
                
        hdu.header.append("BCKGTRANS")
        if data is not None:
            hdu.header['OFFAXISTRANS']=thruput*mask_trans
            hdu.header['BCKGTRANS']=thruput*mask_trans
        else:
            hdu.header['OFFAXISTRANS']=thruput
            hdu.header['BCKGTRANS']=thruput
        if 'APP' in mode:
            hdu.header['OFFAXISTRANS']=hdu.header['OFFAXISTRANS']*app_single_psf # the invididual APP PSF contains a fraction (~48%) of the light passing through the coronagraph optics, but the background still 100% of the light passing through coronagraph.
            psf_OFF=np.median(psf_ON,0)
        hdu.header['CDELT1']=conf['band_specs'][hdu.header['BAND']]['pscale']*0.000277778/1000.
        hdu.header['CDELT2']=conf['band_specs'][hdu.header['BAND']]['pscale']*0.000277778/1000. # mas to degrees
        print(conf)
        hdu.header['EXPTIME']=dit #conf['dit']
        hdu.header["CTYPE1"]  = 'RA---TAN'
        hdu.header["CTYPE2"]  = 'DEC--TAN'
        hdu.header['CUNIT1']='deg'
        hdu.header['CUNIT2']='deg'
        if conf['ScopeSim_LMS']==False:
            hdu.header['BUNIT']="photons/s" # BUNIT overridden by use of Vega spectrum in scopesim
        print("band is ",hdu.header['BAND'])
        if hdu.header['BAND'] in ["N1","N2","N1a","N2a"]:
            hdu.header['CRPIX1']=163
            hdu.header['CRPIX2']=163
        else:
            hdu.header["CRPIX1"]=202
            hdu.header["CRPIX2"]=202
        hdu.header["CRVAL1"]=0
        hdu.header["CRVAL2"]=0
        
        if conf['ScopeSim_LMS']==True:
            hdu.header['NAXIS1']=403#2048
            hdu.header['NAXIS2']=403#2048##conf['band_specs'][hdu.header['BAND']]['pscale']*0.000277778/1000. # mas to degrees
            hdu.header['NAXIS3']=20
            #hdu.header['CDELT1']=5.78*0.000277778/1000.
            #hdu.header['CDELT2']=5.78*0.000277778/1000.##conf['band_specs'][hdu.header['BAND']]['pscale']*0.000277778/1000. # mas to degrees
#        print(conf)
            hdu.header['CDELT3']=0.00005
            hdu.header['EXPTIME']=6#0.011#dit#conf['dit']
            hdu.header["CTYPE1"]  = 'RA---TAN'
            hdu.header["CTYPE2"]  = 'DEC--TAN'
            hdu.header['CUNIT1']='deg'
            hdu.header['CUNIT2']='deg'
            hdu.header['CUNIT3']='um'
            hdu.header['BUNIT']="erg/cm^2/s/angstrom" # BUNIT overridden by use of Vega spectrum in scopesim
            hdu.header["CRPIX1"]=202#1024
            hdu.header["CRPIX2"]=202#1024 # this needs to be half of naxis1 otherwise it will put the star at the wrong center
            hdu.header["CRPIX3"]=0
            hdu.header["CRVAL1"]=0
            hdu.header["CRVAL2"]=0
            hdu.header["CRVAL3"]=hdu.header['LAM']/1e-6 # start value of wavelength range in microns
            psf_ON=np.repeat(psf_ON[:,np.newaxis,:,:],hdu.header['NAXIS3'],axis=1)
        
        header=hdu.header
        print(conf)
        psf_ON, psf_OFF = sim_heeps(psf_ON*10**(-0.4*(mag-0.)), psf_OFF*10**(-0.4*(mag-0.)), header,**conf) # note, mask_trans and offaxis_trans not yet transferred to background
        
    else: # if not passed to ScopeSim
        # rescale PSFs to star signal
        star_signal = dit * flux_star * 10**(-0.4*(mag - mag_ref))
        psf_OFF *= star_signal * mask_trans
        psf_ON *= star_signal * mask_trans
        # add background
        bckg_noise = dit * flux_bckg * thruput * mask_trans
        psf_ON += bckg_noise
        # add photon noise ~ N(0, sqrt(psf))
        np.random.seed(seed)
        psf_ON += np.random.normal(0, np.sqrt(psf_ON))

    if data is None:
        mask_trans=1.

    if (verbose is True) and (call_ScopeSim==False):
        print('   offaxis_trans=%3.4f, mask_trans=%3.4f,'%(thruput, mask_trans))
        print('   mag=%s, dit=%3.3f'%(mag, dit))
        print('   star_signal=%3.2E, bckg_noise=%3.2E'%(star_signal, bckg_noise))
    elif (verbose is True) and (call_ScopeSim==True) and (conf['ScopeSim_LMS']==False):
        print('   SCOPESIM MODE: some values biased')
        print('   offaxis_trans=%3.4f, mask_trans=%3.4f,'%(thruput, mask_trans))
        print('   mag=%s, dit=%3.3f'%(mag, dit))
        print('   star_signal=%3.4f, bckg_noise=%3.4f'%(np.sum(psf_OFF)/thruput/mask_trans, np.median(psf_ON)))
    elif (verbose is True) and (call_ScopeSim==True) and (conf['ScopeSim_LMS']==True):
        print('   SCOPESIM LMS MODE: some values biased')
        print('   offaxis_trans=%3.4f, mask_trans=%3.4f,'%(thruput, mask_trans))
        print('   mag=%s, dit=%3.3f'%(mag, dit))
        print('   star_signal=%3.4f, bckg_noise=%3.4f'%(np.sum(psf_OFF)/thruput/mask_trans, np.median(psf_ON)))

    return psf_ON, psf_OFF
