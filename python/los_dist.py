import numpy as np
from astropy import cosmology
import matplotlib.pyplot as plt

zl = 0.86
zs= 1.28

Re_arc = 1.36 # arcsec
Rf_arc = 0.7843
Rf_arr = np.zeros(100)
Rf_arr.fill(Rf_arc)

#Re = np.radians(Re_arc/3600.)*cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h #Mpc


## foreground
z_list = np.linspace(0,zl,100)
## background
z2_list = np.linspace(zl,zs,100)

dist_ratio = cosmology.Planck13.comoving_distance(z_list)/cosmology.Planck13.comoving_distance(zl)
Dls = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(zl))
Dxs = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(z2_list))
dist2_ratio = (Dxs/Dls)

plt.plot(z_list,Rf_arr,color='b',label='observed position')
plt.plot(z2_list,Rf_arr,color='b')
plt.plot(z_list,Re_arc*dist_ratio,color='r',label='los position')
plt.plot(z2_list,Re_arc*dist2_ratio,color='r')
plt.legend(loc=2)
plt.savefig('B2045_los_test.png')