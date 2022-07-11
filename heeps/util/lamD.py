import astropy.units as u

def get_lamD_pix(lam=3.81e-6, diam_ext=36.905, diam_nominal=38.542,
        ls_dRext=0.0209, pscale=5.47, **conf):
    
    # 1 lambda/D in mas
    lamD = get_lamD_mas(lam=lam, diam_ext=diam_ext, diam_nominal=diam_nominal, 
        ls_dRext=ls_dRext)
    # 1 lambda/D in pixels
    lamD /= pscale
    
    return lamD

def get_lamD_mas(lam=3.81e-6, diam_ext=36.905, diam_nominal=38.542,
        ls_dRext=0.0209, **conf):
    
    # 1 lambda/D in mas
    diamLS = diam_ext - diam_nominal*ls_dRext
    lamD = lam/diamLS*u.rad.to('mas')

    return lamD