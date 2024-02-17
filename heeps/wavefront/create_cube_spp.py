from heeps.util.multiCPU import multiCPU

band = 'L'
mode = 'SPP'
ndet = 293
cw = 19
cmin = int((ndet-cw)/2)
cmax = int((ndet+cw)/2)

filename = 'onaxis_PSF_%s_%s'%(band, mode) + '%s.fits'
spp_stripe = fits.getdata(filename%'_stripe')

def build_SPP(cmin, cmax, frame):
    frame[cmin:cmax,:] = np.transpose(frame[:,cmin:cmax], (1, 0))
    return frame

case = 'create %s band %s cube'%(band, mode)
spp = multiCPU(build_SPP, case=case,
                    posargs=[cmin, cmax],
                    posvars=[spp_stripe])
print('saving ' + filename%'')
fits.writeto(filename%'', spp, overwrite=True)