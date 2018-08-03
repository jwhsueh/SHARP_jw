import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

'''
This code calculate the half-mode mass of WDM
'''
wdm_m = 2.0 # keV
Owdm = cosmo.Om0
h=cosmo.h
mu=1.12
rho = cosmo.critical_density(0) # g/cm^3
rho = rho.to(u.Msun/u.Mpc**3)

alpha = 0.049*(wdm_m)**(-1.11)*(Owdm/0.25)**0.11*(h/0.7)**1.22*(u.Mpc/u.h) # Mpc h-1
Mfs = 4.0*np.pi/3.0*rho*(alpha/2.0)**3
Mhm = 2.7e3*Mfs

print Mhm*h**3
