import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


hdul = fits.open('../data/SIE_kappa.fits')
kappa_fits = hdul[0].data
#print kappa.shape
kappa = np.array(kappa_fits[500,501:])
#print kappa

re = np.linspace(0.01,5.0,499)
print re[50]

kappa_cir = np.zeros(499)
for i in range(499):
	kappa_a = kappa[i]*np.pi*(re[i]**2-(re[i]-0.01)**2)
	kappa_cir[i] = kappa_a
	if i>0:
		kappa_cir[i] = kappa_cir[i]+kappa_cir[i-1]

	kappa_cir[i] = kappa_cir[i]/(np.pi*re[i]**2)

plt.plot(re[150:],kappa_cir[150:])
plt.show()