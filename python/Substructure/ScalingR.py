"""This code compute scaling radius of a NFW halo from velocity dispersion"""

import numpy as np
import AngularDistance as distance

## parameters
q=0.3
b=0.5 # arcsec

zl=0.4
zs=1.3

c = 3e8 # m/s
G = 6.67e-11 # m^3/kg/s^2


## cosmology parameter (cosmopara class)
class cosmopara:
	OM = 0.27
	h = 0.68

## velocity dispersion

Ds = distance.angular_distance(cosmopara, zs)
Dl = distance.angular_distance(cosmopara, zl)
Dls = Ds - Dl

b = np.radians(b/3600)

sig = b * np.sqrt(q)/(4*np.pi)*(Ds/Dls)*c**2
sig = np.sqrt(sig)

## M200, c200 & scale length

Hz = 100*cosmopara.h*distance.Ez(cosmopara,zl) # km/s/Mpc
Hz = Hz*1e3/3.08e22 # s-1
#print Hz
#Hz = distance.Ez(cosmopara,zl) # h

#M_200 = sig**3*(G/100/Hz**2)**(0.5)*(2/G)**(3/2) # kg
M_200 = sig**3*(G/Hz**2)**(0.5)*(2/G)**(3./2.) # kg
M_200_sun = M_200/2e30 # M_sun

print M_200_sun

c_200 = 5.71*(M_200_sun/2e12*cosmopara.h)**(-0.084)*(1+zl)**(-0.47) # no unit

print c_200

rho_c = 3*Hz**2/(8*np.pi*G) # kg/m^3

rs = (3*M_200/8/np.pi/rho_c/c_200**3)**(1.0/3.0) # m * h^-1
rs = rs/3.08e22/cosmopara.h/100 # convert to Mpc
print rs
rs = np.degrees(rs/Dl)*3600 # arcsec 

print rs
