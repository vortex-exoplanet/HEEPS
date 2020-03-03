import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from astropy.io import fits
import astropy.units as u


# HEEPS inputs
case_name = 'drift' #all_effects_2048 #scao_only_2048
mode = 'RAVC'
scao = 'compass3600s_samp300ms'    #'compass600s_samp100ms'  #'compass3600s_samp300ms'
band = 'M'
band_bckg = '%s_bckg0'%band
#band_bckg = '%s_bckg1'%band #'L_bckg1'
dec = '-2.47'
diam = 37                       # diameter in m
r_CLC = 2.5                     # radius of CLC occulter in lam/D


# general figure parameters
figsize = (12, 4)#None#
figsize_log = (12, 6)
markevery = 0.12
markersize = 4
title_pad = -20
title_box = dict(edgecolor='w', facecolor='w')#, alpha=0.7)
loglog = True
add_top_xticks = True
xticks_log = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
ylim = (1e-8,1e-3)#(1e-8,1e-1)

# specific figure parameters:
# - single plot: '1' or '1b' (2 columns, e.g. with/without shot noise)
# - double plot: '2' or '2b' (2 columns, e.g. with/without shot noise)
figtypes = ['1']
figtype = figtypes[0]
hspace = 0.02
label_b = 'without shot noise'
xlim_log = (1, 30)#(1, 50)#
ylims = [(1e-6,1e-0), (1e-8,9e-4)]
colors = ['k', 'C1', 'C0', 'C3', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'r', 'v', 'b']
markers = ['+', 'o', 'x', 'd', 'v', '^', '>', '<', '*', 's', 'P', 'X']
linestyles = ['-']*14
two_cols = False
two_cols_arr = np.array([0,2,4,6,8,1,3,5,7,9])
linestyle_b = '--'
linewidth_b = 0.8
legend_size = 'medium'
if case_name is 'static_ncpa_noLSF':
    colors = ['k', 'C1', 'C1', 'C0', 'C0', 'C2', 'C2']
    markers = ['+', 'o', 'o', 'x', 'x', 'd', 'd']
    linestyles = ['-', '-', '--', '-', '--', '-', '--']
elif case_name is 'seeing':
    scao = 'compass600s_samp100ms'
    ylims = [(1e-7,1e-3),(1e-8,9e-5)]
#    two_cols = True
else:
    ylims = [ylim] if figtype in ['1', '1b'] else [(1e-6,1e-2),(1e-8,9e-4)]
    legend_size = 'small'


# choose figure types:
#two_cols_arr = np.array([0,2,4,1,3,5])


# general filename
filename = 'cc_%s_dec%s_%s_%s_%s'\
        %(scao, dec, band_bckg, mode, case_name) + '_%s.fits'
# filenames: e.g. 1=adi, 2=raw contrast
filename1 = filename%('ADI' if figtype in ['1', '1b'] else 'raw')
filename2 = filename%'ADI'
# secondary filenames (e.g. with/without shot noise)
filename1b = None
filename2b = None
# saved figure file names
suffix = '_log' #'_log' if loglog is True else '_lin'
savename = filename1[:-5] + '_%s' + '%s.png'%suffix


# load filename1 data 
data1 = fits.getdata(filename1)
HDR1 = fits.open(filename1)[0].header
# assume same lambda/D for all curves
lam = HDR1['lam']
lamD = (lam/diam*u.rad).to('arcsec').value
# load other data tables, if present
data2 = fits.getdata(filename2) if filename2 else None
HDR2 = fits.open(filename2)[0].header if filename2 else None
data1b = fits.getdata(filename1b) if filename1b else None
data2b = fits.getdata(filename2b) if filename2b else None


# loop through figure types
for figtype in figtypes:
    
    # single plot
    if figtype in ['1', '1b']:
        fig = plt.figure(figsize=figsize)
        axes = [fig.gca()]
        axes[0].set_title(HDR1['title'], pad=title_pad, bbox=title_box)
        HDRs = [HDR1]
        dataAs = [data1]
        dataBs = [data1b]
    
    # double plot
    elif figtype in ['2', '2b']:
        fig = plt.figure(figsize=figsize_log)
        fig.subplots(2, 1, sharex=True)
        fig.subplots_adjust(hspace=hspace)
        axes = fig.axes
        axes[0].set_title(HDR1['title'], pad=title_pad, bbox=title_box)
        axes[1].set_title(HDR2['title'], pad=title_pad, bbox=title_box)
        HDRs = [HDR1, HDR2]
        dataAs = [data1, data2]
        dataBs = [data1b, data2b]
    
    # draw curves (single or double plot)
    for ax, dataA, dataB, HDR, ylim in zip(axes, dataAs, dataBs, HDRs, ylims):
        # initialize handle of curves (for legend)
        h = []
        # get x-values in lam/D
        xA = dataA[0]
        # loop through y-values
        for i, yA in enumerate(dataA[1:]):
            # get label
            label = HDR['label' + str(i+1)]
            # if CLC mode, create mask to clip data to CLC occulter radius
            r = r_CLC if 'CLC' in label else 0
            mask = (xA > r)
            h.append(ax.plot(xA[mask], yA[mask], label=label, \
                    color=colors[i], marker=markers[i], \
                    markersize=markersize, markevery=markevery, \
                    linestyle=linestyles[i])[0])
            # add secondary curves
            if 'b' in figtype:
                yB = dataB[i+1]
                h.append(ax.plot(xA[mask], yB[mask], label=label_b, \
                        color=colors[i], marker=markers[i], \
                        markersize=markersize, markevery=markevery, \
                        linestyle=linestyle_b, linewidth=linewidth_b)[0])
            # set y limit
            ax.set_ylim(ylim)
    
    # customize figure axes
    for j, (ax, HDR) in enumerate(zip(axes, HDRs)):
        if loglog is True:
            ax.loglog()
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.set_xticks([1,2,5,10,20,50])
            ax.set_xlim(xlim_log)
        else:
            ax.set_yscale('log')
            ax.set_xlim(left=0)
        if j == len(axes)-1:
            ax.set_xlabel(HDR['xlabel'])
        ax.set_ylabel(HDR['ylabel'])
        ax.grid()
        ax.grid(which='minor', linestyle=':')
        # add arcsec ticks on top
        if add_top_xticks is True:
            if j == 0:
                ax1_xlim = np.array(ax.get_xlim())
                ax2_xlim = ax1_xlim*lamD
                ax2 = ax.twiny()
                ax2.set_xlim(ax2_xlim)
                if loglog is True:
                    ax2.set_xscale('log')
                    xt2 = [x for x in xticks_log if ax2_xlim[0] <= x <= ax2_xlim[1]]
                    ax2.set_xticks(xt2)
            str_ticks = [str(round(x,2))+'"' for x in ax2.get_xticks()]
            ax2.set_xticklabels(str_ticks)
#            ax2.text(0.01, 1.04, '%s-band %s'%(band, mode), transform=ax2.transAxes)
#        horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

        # sort legend in two columns
        if two_cols is True:
            curves = np.array(h)
#            curves = np.array(h)[two_cols_arr] \
#                if band in ['L','M'] else np.array(h)[np.array([0,2,4,1,3,5])]
            labs = [c.get_label() for c in curves]
#            labs[0] = 'ELT'
            ax.legend(curves, labs, loc='upper right', ncol=2)
        else:
            ax.legend(loc='upper right', prop={'size': legend_size})
    # save figure
    if True:
        plt.show(block=False)
        plt.savefig(savename%figtype, dpi=300, transparent=True)
        #plt.close()
