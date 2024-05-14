from heeps.util.multiCPU import multiCPU
from imutils import rotate

band = 'L'
mode = 'SPP'
ndet = 293
cw = 23#19
cmin = int((ndet-cw)/2)
cmax = int((ndet+cw)/2)
rot = -10

filename = 'onaxis_PSF_%s_%s'%(band, mode) + '%s.fits'
spp_stripe = fits.getdata(filename%'_stripe')
case = 'create %s band %s cube'%(band, mode)

mask = np.zeros((ndet,ndet))
mask[:,cmin:cmax] = 1
mask_rot = rotate(mask, rot)

def build_SPP_rot(mask_rot, frame):
    stripe = np.rot90(frame*mask_rot)
    frame *= (1-np.rot90(mask_rot))
    frame += stripe
    return frame
spp = multiCPU(build_SPP_rot, case=case,
                    posargs=[mask_rot],
                    posvars=[spp_stripe])

#def build_SPP(cmin, cmax, frame):
#    frame[cmin:cmax,:] = np.transpose(frame[:,cmin:cmax], (1, 0))
#    return frame
#spp = multiCPU(build_SPP, case=case,
#                    posargs=[cmin, cmax],
#                    posvars=[spp_stripe])

print('saving ' + filename%'')
fits.writeto(filename%'', spp, overwrite=True)