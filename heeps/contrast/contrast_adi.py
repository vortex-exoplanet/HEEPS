import matplotlib.pyplot as plt
from astropy.io import fits

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
bands = ['L', 'M', 'N1', 'N2']#['N1', 'N2']# 
mag = 5# -1.6#
bckg = False

# modes per band
modes = {'L': ['ELT', 'CVC', 'RAVC'],
         'M': ['ELT', 'CVC', 'RAVC'],
         'N1': ['ELT', 'CVC'],
         'N2': ['ELT', 'CVC']}
colors = ['k', 'C0', 'C1']

# file names
loadname = 'cc_compass600s_samp100ms_ADI3600s_samp100ms_avg0ms_dec-2.47deg_%s_mag%s_bckg%s_%s.fits'
savename = 'cc_adi_%s_mag%s_bckg%s_%s.png'

def savefig(fignum, band, suffix, mag=mag, bckg=bckg, savename=savename):
    plt.figure(fignum)
    plt.yscale('log')
    plt.grid()
    plt.grid(which='minor', linestyle=':')
    plt.xlabel("Angular separation [arcsec]")
    plt.ylabel(r"5-$\sigma$ sensitivity")
    wwo = {1:'with', 0:'without'}
    plt.title(r"Star mag %s = %s, %s background"%(band, mag, wwo[bckg]))
    plt.legend()
    plt.xlim(left=0)
    plt.ylim(1e-8, 1e-2)
    plt.show(block=False)
    plt.savefig(savename%(band, mag, int(bckg), suffix), dpi=300, transparent=True)
    plt.close()

for band in bands:
    for i, mode in enumerate(modes[band]):
        data = fits.getdata(loadname%(band, mag, int(bckg), mode))
        # normal Gaussian distribution
        plt.figure(1)
        plt.plot(data[:,4],data[:,0], color=colors[i], label=mode)
        # Student's distribution
        plt.figure(2)
        plt.plot(data[:,4],data[:,1], color=colors[i], label=mode)
    
    # save figures
    savefig(1, band, 'normal')
    savefig(2, band, 'student')
