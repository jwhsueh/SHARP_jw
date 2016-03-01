""" This code test different profiles in 2D projected position """

import numpy as np
import NFWprofile as NFW
import scipy.stats
import scipy.integrate

import matplotlib.pyplot as plt

## lenspara

class lenspara:
	q=0.3
	b=0.5 # arcsec

	zl=0.4
	zs=1.3

	xc = 0.0 # centroid position
	yc = 0.0

##------------
##-----field of view--------

## 2x Einstein radius (b)

x_lim, y_lim = lenspara.xc+np.array([-lenspara.b,lenspara.b]), lenspara.yc+np.array([-lenspara.b,lenspara.b])



## ------NFW profile para setting------##

rs = NFW.scaleR(lenspara)
r_end = 10*rs

## -------------

n_draw = 100000

rand = np.random.rand(n_draw)

r_d1 = NFW.inverse_cdf(rand,rs,r_end)
r_d2 = rand*r_end # uniform distribution

## 3D->2D projection

# draw angle theta[0~pi] & phi[0~2*pi]
# phi is on the plane perpendicular to our line of sight

theta_d = rand*np.pi
phi_d = rand*2*np.pi

# x,y coordinate of each pt

## 3D profile

#x_d1 = lenspara.xc+r_d1*np.sin(theta_d)*np.cos(phi_d)
#y_d1 = lenspara.yc+r_d1*np.sin(theta_d)*np.sin(phi_d)

rp_d1 = r_d1*np.sin(theta_d)

plt.hist(rp_d1, bins=np.linspace(0,2,200))

#plt.scatter(x_d1,y_d1)
plt.xlim(0,2*lenspara.b)
plt.xlabel('distance to center in arc sec (b = 0.5")')
plt.ylabel('number counts')
plt.title('NFW profile 2D distribution (N_tot = 100000)')
#plt.ylim(y_lim)
plt.savefig('NFW_2D.png')
#plt.show()