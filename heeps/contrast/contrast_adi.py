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
save2hdu = False
two_cols = True
add_top_xticks = True
xticks_log = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
ylim1 = (1e-7, 1e-2)
ylim2 = (3e-8, 3e-3)#(1e-8, 1e-2)#(2e-7, 3e-3)#
xlabel = r"Angular separation $[\lambda/D]$"
ylabel = r"5-$\sigma$ sensitivity (star mag %s = %s)"

#select figure types: 0=without bckg; 1=with bckg; 2=with and without; 3=double plot
figtypes = [2]#range(4) #all
# add optional suffix
suffix = '_test'
suffix = '_log' if loglog is True else '_lin'

# modes per band
colors = ['k', 'C1', 'C0', 'C3', 'C2']
markers = ['+', 'o', 'x', 'd', 'v']
band_specs = {'L': {'lam': 3.8e-6,
                    'mag': 6,
                  'modes': ['ELT', 'RAVC', 'CVC', 'APP', 'CLC'],
                 'colors': ['k', 'C2', 'C1', 'C3', 'C0'],
                'markers': ['+', 'v', 'o', 'd', 'x']},
              'M': {'lam': 4.8e-6,
                    'mag': 6,
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
savename_fits = 'cc_adi_%s_mag%s_%s_bckg%s' #+ '_%s'%(distrib)

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
            axes[0].set_title(r"Star mag %s=%s with shot noise"%(band, mag), pad=title_pad)
            axes[1].set_title(r"%s band, no shot noise"%(band), pad=title_pad)
#            axes[0].set_ylim(ylim1)
#            axes[1].set_ylim(ylim2)
            bckgs = [True, False]
            linestyles = ['-', '-']
            linewidths = [1, 1]
            labels = ['%s','%s']
        else:
#            plt.ylim(ylim2)
            wwo = {True:'with', False:'without'} # labels: with or without shot noise
            if figtype in [0,1]:
                axes = [fig.gca()]
                bckgs = [bool(figtype)]
                linestyles = ['-']
                linewidths = [None]
                labels = ['%s'+' %s shot noise'%wwo[bool(figtype)]]
                if figtype is 0:
                    plt.title(r"%s band, no shot noise"%(band))
                else:
                    plt.title(r"Star mag %s = %s"%(band, mag))
            else:
                axes = [fig.gca(), fig.gca()]
                bckgs = [True, False]
                linestyles = ['-', '--']
                linewidths = [None, 0.8]
                labels = ['ELT + %s', 'without shot noise']
        # initialize handle of curves (for legend)
        h = []
        # loop modes for each band figure
        for i, (mode, color, marker) in enumerate(zip(modes, colors, markers)):
            # plot 1 =  with background, plot 2 = without background
            for ax, bckg, linestyle, linewidth, label in \
                    zip(axes, bckgs, linestyles, linewidths, labels):
                data_mag = 5 if band in ['L', 'M'] and bckg is False else mag
                data = fits.getdata(loadname%(band,data_mag, int(bckg), mode))
                x = data[:,4]/(lam/diam*u.rad).to('arcsec').value
                y = data[:,1] if student is True else data[:,0]
                # clip CLC to occulter radius
                if mode == 'CLC':
                    mask = np.where(x>r_CLC)[0]
                    x = x[mask]
                    y = y[mask]
                # add mode to label, when necessary
                try:
                    label = label%mode
                except:
                    pass
                # create plots
                h.append(ax.plot(x, y, linestyle=linestyle, linewidth=linewidth, \
                        label=label, color=color, marker=marker, \
                        markersize=markersize, markevery=markevery)[0])
                # save (x,y) for each mode in fits.hdu (optional)
                if save2hdu is True:
                    date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                    header = fits.Header({'xlabel':xlabel,'ylabel':ylabel%(band,mag),\
                            'date':date, 'wavelnth':lam, 'band':band, \
                            'mag':data_mag, 'bckg':bckg, 'mode':mode})
                    hdu = fits.PrimaryHDU((x,y), header=header)
                    hdu.writeto(savename_fits%(band, data_mag, mode, int(bckg)) \
                            + '.fits', overwrite=True)
        # customize figure axes (keep only 1 axis for figtype 2)
        if figtype is 2:
            axes = [axes[0]]
        for i,ax in enumerate(axes):
            if loglog is True:
                ax.set_xlim(left=1)
                ax.loglog()
                ax.xaxis.set_major_formatter(ScalarFormatter())
                ax.set_xticks([1,2,5,10,20,50])
            else:
                ax.set_xlim(left=0)
                ax.set_yscale('log')
            if i == len(axes)-1:
                ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel%(band,mag))
            ax.grid()
            ax.grid(which='minor', linestyle=':')
            # add arcsec ticks on top
            if add_top_xticks is True:
                ax1_xlim = np.array(ax.get_xlim())
                ax2_xlim = ax1_xlim*(lam/diam*u.rad).to('arcsec').value
                ax2 = ax.twiny()
                ax2.set_xlim(ax2_xlim)
                if loglog is True:
                    ax2.set_xscale('log')
                    xt2 = [x for x in xticks_log if ax2_xlim[0] <= x <= ax2_xlim[1]]
                    ax2.set_xticks(xt2)
                str_ticks = [str(round(x,2))+'"' for x in ax2.get_xticks()]
                ax2.set_xticklabels(str_ticks)
            # sort legend in two columns
            if two_cols is True:
                curves = np.array(h)[np.array([0,2,4,6,8,1,3,5,7,9])] \
                    if band in ['L','M'] else np.array(h)[np.array([0,2,4,1,3,5])]
                labs = [c.get_label() for c in curves]
                labs[0] = 'ELT'
                ax.legend(curves, labs, loc='upper right', ncol=2)
            else:
                ax.legend(loc='upper right')
        # save figure
        if True:
            plt.show(block=False)
            plt.savefig(savename_png%(band, mag, figtype) + '.png', dpi=300, transparent=True)
            plt.close()
