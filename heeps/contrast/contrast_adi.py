import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import time

"""
VIP cc data:
0: sensitivity_gaussian
1: sensitivity_student
2: throughput
3: distance
4: distance_arcsec
5: noise
6: sigma corr
"""

# inputs
bands = ['L', 'M', 'N1', 'N2']
r_CLC = 2.5 # lam/D of CLC occulter
diam = 37
figsize = (12,4)#None#
title_pad = -20
hspace = 0.05
markevery = 0.12
markersize = 4
student = True
loglog = False
ylim1 = (1e-7, 1e-2)
ylim2 = (1e-8, 1e-2)
xlabel = "Angular separation [arcsec]"
ylabel = r"5-$\sigma$ sensitivity"

#select figure types: 0=without bckg; 1=with bckg; 2=with and without; 3=double plot
figtypes = [2]#range(4) #all
# add optional suffix
suffix = '_long'

# modes per band
colors = ['k', 'C1', 'C0', 'C3', 'C2']
markers = ['+', 'o', 'x', 'd', 'v']
band_specs = {'L': {'lam': 3.8e-6,
                    'mag': 5,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'M': {'lam': 4.8e-6,
                    'mag': 5,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'N1':{'lam': 8.7e-6,
                    'mag': -1.6,
                  'modes': ['ELT', 'CVC', 'CLC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']},
              'N2':{'lam': 11.5e-6,
                    'mag': -1.6,
                  'modes': ['ELT', 'CVC', 'CLC'],
                 'colors': ['k', 'C1', 'C0'],
                'markers': ['+', 'o', 'x']}}

# file names
loadname = 'cc_compass3600s_samp300ms_ADI3600s_samp300ms_avg0ms_dec-2.47deg_%s_mag%s_bckg%s_%s.fits'
distrib = 'student' if student is True else 'normal'
savename_png = 'cc_adi_%s_mag%s_figtype%s' + '_%s%s'%(distrib, suffix)
savename_fits = 'cc_adi_%s_mag%s_%s' #+ '_%s'%(distrib)

# loop figure types
for figtype in figtypes:
    # loop bands and create one figure per band
    for band in bands:
        # get useful specs
        lam = band_specs[band]['lam']
        mag = band_specs[band]['mag']
        modes = band_specs[band]['modes']
        colors = band_specs[band]['colors']
        markers = band_specs[band]['markers']
        # create figure type: 0=without bckg; 1=with bckg; 2=with and without; 3=double plot
        fig = plt.figure(figsize=figsize)
        if figtype is 3:
            fig.subplots(2, 1, sharex=True)
            fig.subplots_adjust(hspace=hspace)
            axes = fig.axes
            axes[0].set_title(r"Star mag %s=%s with background"%(band, mag), pad=title_pad)
            axes[1].set_title(r"%s band, no background"%(band), pad=title_pad)
#            axes[0].set_ylim(ylim1)
#            axes[1].set_ylim(ylim2)
            bckgs = [True, False]
            linestyles = ['-', '-']
            linewidths = [1, 1]
            labels = ['%s','%s']
        else:
            plt.title(r"Star mag %s = %s"%(band, mag))
#            plt.ylim(ylim2)
            wwo = {False:'without', True:'with'} # labels: with or without background
            if figtype in [0,1]:
                axes = [fig.gca()]
                bckgs = [bool(figtype)]
                linestyles = ['-']
                linewidths = [None]
                labels = ['%s'+' %s background'%wwo[bool(figtype)]]
                if figtype is 0:
                    plt.title(r"%s band, no background"%(band))
            else:
                axes = [fig.gca(), fig.gca()]
                bckgs = [True, False]
                linestyles = ['--', '-']
                linewidths = [0.8, None]
                labels = ['%s with background', '%s without background']
        # loop modes for each band figure
        for i, (mode, color, marker) in enumerate(zip(modes, colors, markers)):
            # plot 1 =  with background, plot 2 = without background
            for ax, bckg, linestyle, linewidth, label in \
                    zip(axes, bckgs, linestyles, linewidths, labels):
                data = fits.getdata(loadname%(band, mag, int(bckg), mode))
                x = data[:,4]
                y = data[:,1] if student is True else data[:,0]
                # clip CLC to occulter radius
                if mode == 'CLC':
                    r_mask = (r_CLC*lam/diam*u.rad).to('arcsec').value
                    mask = np.where(x>r_mask)[0]
                    x = x[mask]
                    y = y[mask]
                ax.plot(x, y, linestyle=linestyle, linewidth=linewidth, \
                        label=label%mode, color=color, marker=marker, \
                        markersize=markersize, markevery=markevery)
                # save fits.hdu for each mode
                date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                header = fits.Header({'xlabel':xlabel,'ylabel':ylabel,\
                        'date':date, 'wavelnth':lam, 'band':band, 'mag':mag, 'mode':mode})
                hdu = fits.PrimaryHDU((x,y), header=header)
                hdu.writeto(savename_fits%(band, mag, mode) + '.fits', overwrite=True)
        # customize figure axes (keep only 1 axis for figtype 2)
        if figtype is 2:
            axes = [axes[0]]
        for i,ax in enumerate(axes):
            if loglog is True:
                ax.loglog()
                ax.xaxis.set_major_formatter(ScalarFormatter())
            else:
                ax.set_yscale('log')
                ax.set_xlim(left=0)
            if i == len(axes)-1:
                ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid()
            ax.grid(which='minor', linestyle=':')
            ax.legend()
        # save figure
        if True:
            plt.show(block=False)
            plt.savefig(savename_png%(band, mag, figtype) + '.png', dpi=300, transparent=True)
            plt.close()
