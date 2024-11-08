"""
ScopeSim processing of a HCI (RAVC) data cube produced by HEEPS
The cube contains the temporal evolution of the coronagraphic PSF using
the METIS RAVC in the L band. ScopeSim uses the cube as the input source
and produced detector images with background and detector noise components.
"""
import os
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import matplotlib.pyplot as plt
import scopesim as sim
import scipy.ndimage
from heeps.util.save2fits import save2fits


#sim.download_package(['instruments/METIS', 'telescopes/ELT', 'locations/Armazones'])

from synphot import SourceSpectrum

# Change if you want to use another cube
#FNAME = "onaxis_PSF_L_RAVC.fits"

WORKDIR = "./"
PKGSDIR = "./inst_pkgs/"

# This sets the path to the METIS configuration. Change the path to
# match your setup.
PKGS_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(PKGSDIR)))
sim.rc.__config__['!SIM.file.local_packages_path'] = PKGS_DIR
    

list_lm=[
         "HCI_L_long",
         "HCI_L_short",
         "Lp",
         "short-L",
         "L_spec",
         "Mp",
         "M_spec",
         "Br_alpha",
         "Br_alpha_ref",
         "PAH_3.3",
         "PAH_3.3_ref",
         "CO_1-0_ice",
         "CO_ref",
         "H2O-ice",
         "IB_4.05",
         "HCI_M"]
def sim_heeps(psf_ON, psf_OFF, header, **conf):
    """Do the simulation"""

    offaxistrans=header['OFFAXISTRANS'] # this number transfers the transmission of the coronagraphic optics as seen by an off-axis source, per PSF (in case of APP).  
    bckgtrans=header['BCKGTRANS'] # same but for a background source, the difference is when the APP is used. which only contains 48% of the off-axis flux per PSF
#    print(offaxistrans)
    
    # pre-compensate PSF_OFF and PSF_ON for bckgtrans
    psf_OFF=psf_OFF/bckgtrans 
    psf_ON=psf_ON/bckgtrans
    
    
    # Exposure parameters
    dit = header['EXPTIME']
#    header['FILTER']="HCI_L_short"
    #header['FILTER']='PAH_3.3' OVERRIDE FOR TESTING
    ndit = 1
    # Set up the instrument
    #cmd = sim.UserCommands(use_instrument='METIS', set_modes=['img_lm'])
    #print(conf)
    if conf['ScopeSim_LMS']==False: # for IMG_LM or IMG_N modes
        naxis1, naxis2, naxis3 = (header['NAXIS1'],
                              header['NAXIS2'],
                              header['NAXIS3'])
        if header['FILTER'] in list_lm: # IMG_LM mode
            cmd=sim.UserCommands(use_instrument="METIS", set_modes=["img_lm"],properties={"!OBS.filter_name": header['FILTER'],"!OBS.auto_exposure.fill_frac":1e6,"!OBS.exptime":dit,"!OBS.ndit":ndit,"!OBS.dit": dit,"!OBS.detector_readout_mode":"fast","!OBS.pupil_transmission":bckgtrans})
            cmd2=sim.UserCommands(use_instrument="METIS", set_modes=["img_lm"],properties={"!OBS.filter_name": header['FILTER'],"!OBS.auto_exposure.fill_frac":1e20,"!OBS.exptime":dit,"!OBS.ndit":ndit,"!OBS.dit": dit,"!OBS.detector_readout_mode":"fast","!OBS.pupil_transmission":bckgtrans})
        else: # this is IMG_N
            cmd=sim.UserCommands(use_instrument="METIS", set_modes=["img_n"],properties={"!OBS.filter_name": header['FILTER'],"!OBS.auto_exposure.fill_frac":1e6,"!OBS.exptime":dit,"!OBS.ndit":ndit,"!OBS.dit": dit,"!OBS.detector_readout_mode":"high_capacity","!OBS.pupil_transmission":bckgtrans})
            cmd2=sim.UserCommands(use_instrument="METIS", set_modes=["img_n"],properties={"!OBS.filter_name": header['FILTER'],"!OBS.auto_exposure.fill_frac":1e20,"!OBS.exptime":dit,"!OBS.ndit":ndit,"!OBS.dit": dit,"!OBS.detector_readout_mode":"high_capacity","!OBS.pupil_transmission":bckgtrans})
    elif conf['ScopeSim_LMS']==True:
        naxis1, naxis2, naxis3,naxis4 = (header['NAXIS1'],
                              header['NAXIS2'],
                              header['NAXIS3'],np.shape(psf_ON)[0])
        if header['FILTER'] in ["HCI_L_long","CO_ref"]:
            cmd = sim.UserCommands(use_instrument='METIS', set_modes=['lms'],properties={"!OBS.exptime": dit,
                                 "!SIM.spectral.spectral_resolution": 2000,
                                 "!SIM.spectral.spectral_bin_width": 1e-3,'!OBS.wavelen':header['lam']/1e-6,"!OBS.pupil_transmission":bckgtrans,"!OBS.auto_exposure.fill_frac":1e6})
            #cmd2=sim.UserCommands(use_instrument="METIS", set_modes=["img_lm"],properties={"!OBS.filter_name": header['FILTER'],"!OBS.auto_exposure.fill_frac":1e20,"!OBS.exptime":dit,"!OBS.ndit":ndit,"!OBS.dit": dit,"!OBS.detector_readout_mode":"fast"})
    # somewhere we need to set the cold stop pupil to the off-axis transmission to reduce the background flux correctly
    # similarly we need to rescale the on-axis and off-axis PSFs to pre-correct for that factor.
    
    # Set up detector size to match the input cube - using the full
    # METIS detector of 2048 x 2048 pixels would require too much memory
    
    if conf['ScopeSim_LMS']==False:
        cmd["!DET.width"] = naxis1
        cmd["!DET.height"] = naxis2

    metis_on = sim.OpticalTrain(cmd)
    metis_on['psf'].include=False
    if conf['ScopeSim_LMS']==True: # LMS mode
        
        nplanes = 12000
    #nplanes = naxis3  # for full cube
        #naxis4=np.shape(psf_ON)[0]
        print(psf_ON.shape,"psf shape")
        nplanes = min(nplanes, naxis4)
        #psf_ON_out=np.zeros([nplanes,naxis3,naxis1,naxis2])
        for n in np.arange(nplanes):
            new_wl=header['CRVAL3']+np.arange(naxis3)*header['CDELT3']
            print(new_wl)
            vega=fits.getdata(WORKDIR+"alpha_lyr_stis_008.fits")
            #mag=5
            vega_flux=vega["FLUX"] # angstrom to micron
            vega_wave=vega['WAVELENGTH']/10000. # ergs/cm^2/s/angstrom

            f_flux=scipy.interpolate.interp1d(vega_wave,vega_flux)
            new_flux=f_flux(new_wl)
            #new_flux[500:600]=0
            hdu=fits.PrimaryHDU(data=new_flux[:,np.newaxis,np.newaxis]*psf_ON[n,:,:,:],header=header)
            src_on=sim.source.source.Source(cube=hdu,flux=0*u.mag)
            metis_on.observe(src_on,update=True)
            result = metis_on.readout(detector_readout_mode="auto")[0]
            fig,ax=plt.subplots(2,2)
            print(result[1].header)
            ax[0,0].imshow(result[2].data,origin="lower")

            ax[1,0].imshow(result[3].data,origin="lower")

            ax[0,1].imshow(result[1].data,origin="lower")

            ax[1,1].imshow(result[4].data,origin="lower")
            plt.show()
            plt.close()
        
    elif conf['ScopeSim_LMS']==False:
        metis_sum = sim.OpticalTrain(cmd)
        metis_sum['psf'].include=False    
        metis_off = sim.OpticalTrain(cmd2)
        metis_off['psf'].include=False
        metis_off['readout_noise'].include=False
        metis_off['shot_noise'].include=False
        metis_bkg = sim.OpticalTrain(cmd2)
        metis_bkg['psf'].include=False
        metis_bkg['readout_noise'].include=False
        metis_bkg['shot_noise'].include=False
        #print(metis_bkg.effects)
    # We attach the spectrum of Vega to the image
        spec = SourceSpectrum.from_file(WORKDIR+"alpha_lyr_stis_008.fits")

    # Loop over the input cube
        nplanes = 12000
    #nplanes = naxis3  # for full cube
        nplanes = min(nplanes, naxis3)
        #nplanes=4
        psf_ON_out=np.zeros([nplanes,naxis1,naxis2])
        psf_SUM_out=np.zeros([nplanes,naxis1,naxis2])
        psf_SUM=np.zeros([nplanes,naxis1,naxis2])
        psf_OFF_out=np.zeros([naxis1,naxis2])
    
        hourangles=np.linspace(-0.5/24.*360.,0.5/24.*360.,nplanes)
    
        lat = -24.59            # deg
        dec = -51.0665168055                # deg -51 == beta pic
        hr = np.deg2rad(hourangles)
        dr = np.deg2rad(dec)
        lr = np.deg2rad(lat)
# parallactic angle in deg
        pa = -np.rad2deg(np.arctan2(-np.sin(hr), np.cos(dr)*np.tan(lr)
     - np.sin(dr)*np.cos(hr)))
        pa = (pa + 360)%360 
        pa_rad = np.deg2rad(pa)
        ang_sep=100.
        xy = np.array([ang_sep*np.cos(np.pi/2.+pa_rad), ang_sep*np.sin(np.pi/2.+pa_rad)]).T
        for n in np.arange(nplanes):
            if header['FILTER'] in list_lm: # LM band
                psf_SUM[n,:,:]=psf_ON[n,:,:]+scipy.ndimage.shift((10**(-0.4*7.7)*psf_OFF[:,:]),(xy[n,0]/5.47,xy[n,1]/5.47))#+scipy.ndimage.shift((10**(-0.4*7.7)*psf_OFF[:,:]),(xy[n,0]/5.47,xy[n,1]/5.47))
            else: #n-band
                psf_SUM[n,:,:]=psf_ON[n,:,:]+scipy.ndimage.shift((10**(-0.4*5)*psf_OFF[:,:]),(xy[n,0]/5.47,xy[n,1]/5.47))#+scipy.ndimage.shift((10**(-0.4*5)*psf_OFF[:,:]),(xy[n,0]/5.47,xy[n,1]/5.47))
    
        for i in range(nplanes):
            print("Plane", i+1, "/", nplanes)
            print(np.sum(psf_ON[i,:,:]))
            imhdu_on=fits.PrimaryHDU(data=psf_ON[i,:,:],header=header)
            imhdu_off=fits.PrimaryHDU(data=psf_OFF/10000.,header=header) # the off-axis reference should not be saturated; hack to prevent saturation
            imhdu_sum=fits.PrimaryHDU(data=psf_SUM[i,:,:],header=header)
            src_on = sim.Source(spectra=[spec], image_hdu=imhdu_on)
            src_off = sim.Source(spectra=[spec], image_hdu=imhdu_off)
            src_sum = sim.Source(spectra=[spec], image_hdu=imhdu_sum)
            src_bkg = sim.source.source_templates.empty_sky()

            metis_on.observe(src_on,update=True)
            metis_off.observe(src_off,update=True)
            metis_bkg.observe(src_bkg,update=True)
            metis_sum.observe(src_sum,update=True)
            #print(metis_bkg['detector_readout_parameters'].list_modes())
            outhdu = metis_on.readout()[0]
            outhdu2 = metis_off.readout()[0]
            outhdu3 = metis_bkg.readout()[0]
            outhdu4 = metis_sum.readout()[0]
            num=np.int32((naxis1-1)/2)
            print(outhdu[0].header)
            print(outhdu[1].header)
            psf_ON_out[i,:,:] = outhdu[1].data.astype(np.float32)[(1024-num):(1024+num+1),(1024-num):(1024+num+1)]
            psf_OFF_out = outhdu2[1].data[(1024-num):(1024+num+1),(1024-num):(1024+num+1)]-outhdu3[1].data[(1024-num):(1024+num+1),(1024-num):(1024+num+1)]
            psf_OFF_out = psf_OFF_out*10000. # restore flux level to original scale
            psf_SUM_out[i,:,:] = outhdu4[1].data.astype(np.float32)[(1024-num):(1024+num+1),(1024-num):(1024+num+1)]
            newheader_ON=outhdu
            newheader_OFF=outhdu2
            newheader_SUM=outhdu4
    if conf['ScopeSim_LMS']==False:
        newheader_ON[1].data=psf_ON_out
        newheader_OFF[1].data=psf_OFF_out
        newheader_SUM[1].data=psf_SUM_out
        newheader_ON[1].header['CRPIX1']=num+1
        newheader_OFF[1].header['CRPIX1']=num+1
        newheader_SUM[1].header['CRPIX1']=num+1
        newheader_ON[1].header['CRPIX2']=num+1
        newheader_OFF[1].header['CRPIX2']=num+1
        newheader_SUM[1].header['CRPIX2']=num+1

    if True:
        newheader_ON.writeto("scopesim_on_scopesim_{0}_{1}.fits".format(header["MODE"],header["FILTER"]),overwrite=True)
        newheader_OFF.writeto("scopesim_off_scopesim_{0}_{1}.fits".format(header["MODE"],header["FILTER"]),overwrite=True)
        newheader_SUM.writeto("scopesim_sum_scopesim_{0}_{1}.fits".format(header["MODE"],header["FILTER"]),overwrite=True)
    return psf_ON_out, psf_OFF_out

if __name__ == "__main__":
    sim_heeps()
