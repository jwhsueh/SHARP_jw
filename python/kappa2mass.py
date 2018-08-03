import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


hdul = fits.open('../data/B1555_kappa.fits')
kappa_fits = np.array(hdul[0].data)

res = 0.001 #arcsec
d_area = res**2

Re = 1.776548e-01
slen = kappa_fits.shape[0]
cen = slen/2.0
xx,yy =np.meshgrid(np.arange(0,slen),np.arange(0,slen))
xx,yy = xx-cen,yy-cen
rr = np.sqrt(xx**2+yy**2)
mask = rr <= 2.0*Re/res

kappa_tot = np.sum(kappa_fits[mask])*d_area
print kappa_tot
