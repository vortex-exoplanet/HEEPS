#import matplotlib; matplotlib.use('agg') # to run on a headless server
import matplotlib.pyplot as plt
import os.path
from astropy.io import fits
import heeps.util.img_processing as impro
import numpy as np

# inputs
bands = ['Lp', 'Mp', 'N1', 'N2'] #['Lp'] #
pscales = [5, 5, 10, 10] #[5] #
N = 512
(xc,yc) = (int(N/2),int(N/2))
xbin = 1
rim = 250
(xo,yo) = (rim,rim)
ylim = (1e-7, 1e0)

# modes per band
modes = {'Lp': ['VC', 'RAVC'],
         'Mp': ['VC', 'RAVC'],
         'N1': ['VC'],
         'N2': ['VC']}
colors = ['C0', 'C1']
linestyles = ['-','--']

# working repository
#folder = '/mnt/disk4tb/METIS/heeps-analysis/'
folder = '$HOME/INSTRUMENTS/METIS/heeps-analysis/'
path_offaxis = 'offaxis'
path_onaxis = '../cube_COMPASS_20180223_600s_100ms'
path_output = 'output_files'
# absolute paths
folder = os.path.expandvars(folder)
path_offaxis = os.path.join(folder, path_offaxis)
path_onaxis = os.path.join(folder, path_onaxis)
path_output = os.path.join(folder, path_output)
# saved file name
savename = 'cc_raw_%s.png'


def savefig(fignum, band, ylim=ylim, savename=savename):
    plt.figure(fignum)
    plt.yscale('log')
    plt.grid()
    plt.grid(which='minor', linestyle=':')
    plt.xlabel("Angular separation [arcsec]")
    plt.ylabel(r"Raw contrast")
    plt.title(r"%s band"%band)
    plt.legend()
    plt.xlim(left=0)
    plt.ylim(ylim)
    plt.show(block=False)
    plt.savefig(savename%band, dpi=300, transparent=True)
    plt.close()

for band, pscale in zip(bands, pscales):
    x = np.arange(rim+1)*pscale*1e-3*xbin
    plt.figure(1)
    for i, mode in enumerate(modes[band]):
        # off-axis PSF
        psf_OFF = fits.getdata(os.path.join(path_offaxis, 'PSF_%s_%s.fits' \
                %(mode, band)))
        psf_OFF_rim = psf_OFF[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
        # on-axis PSFs (cube)
        psf_ON = fits.getdata(os.path.join(path_onaxis, 'PSF_%s_%s.fits' \
                %(mode, band)))[:15]
        if psf_ON.ndim != 3:
            psf_ON = np.array(psf_ON, ndmin=3)
        # averaged on-axis PSF
        psf_ON_avg = np.mean(psf_ON, 0)
        psf_ON_rim = psf_ON_avg[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
        # radial profiles
        y1 = impro.get_radial_profile(psf_OFF_rim, (xo,yo), xbin)
        y2 = impro.get_radial_profile(psf_ON_rim, (xo,yo), xbin)
        # normalize by the peak of the off-axis PSF
        peak = np.max(y1)
        y1 /= peak
        y2 /= peak
        # clip small/negative values to background level
        bkgd = np.var(psf_ON_avg/peak)
        y1[y1<ylim[0]] = np.min(y1[y1>bkgd])
        y2[y2<ylim[0]] = np.min(y2[y2>bkgd])
        
        # figures
        plt.plot(x, y1, 'k', linestyle=linestyles[i], label='%s off-axis'%mode)
        plt.plot(x, y2, color=colors[i], linestyle=linestyles[i], label='%s on-axis'%mode)
    savefig(1, band)


