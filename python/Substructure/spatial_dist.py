""" This code test different profiles in 2D projected position """

import numpy as np
import matplotlib.pyplot as plt

import NFWprofile as NFW
import MassFunc as MF
import DistanceTool as distance

## lenspara

class lenspara:
	q=0.3
	b=0.5 # arcsec

	zl=0.4
	zs=1.3

	xc = 0.0 # centroid position
	yc = 0.0

## cosmology parameter (cosmopara class)
class cosmopara:
	OM = 0.27
	h = 0.68

##------------
##-----field of view--------

## 2x Einstein radius (b)

x_lim, y_lim = lenspara.xc+np.array([-2.0*lenspara.b,2.0*lenspara.b]), lenspara.yc+np.array([-2.0*lenspara.b,2.0*lenspara.b])

## ------NFW profile para setting------##

rs = NFW.scaleR(lenspara)
r_end = 100*lenspara.b

## -----substructure setting---------##

ml,mu = 1e6, 1e8 # lower & upper limits
f_sub = 0.01 # substructure mass fraction

## -------------

def NFW_2d(n_draw):
#n_draw = 100000

	uni_r = np.random.uniform(n_draw)

	r_d = NFW.inverse_cdf(uni_r,rs,r_end)


## 3D->2D projection

# draw angle theta[0~pi] & phi[0~2*pi]
# phi is on the plane perpendicular to our line of sight

	theta_d = np.random.uniform(0,np.pi,n_draw)
	phi_d = np.random.uniform(0,2.0*np.pi,n_draw)


	#plt.hist(np.degrees(phi_d))

	## 3D to 2D distribution

	rp_d = r_d*np.sin(theta_d) # projected distance to lens center

	x_d = lenspara.xc+rp_d*np.cos(phi_d)
	y_d = lenspara.yc+rp_d*np.sin(phi_d)

	return x_d,y_d,rp_d


def uniform_2d(n_draw):
	x = lenspara.xc+np.random.uniform(x_lim[0],x_lim[1],n_draw)
	y = lenspara.yc+np.random.uniform(y_lim[0],y_lim[1],n_draw)

	r = np.sqrt(x_d**2+y_d**2)

	return x,y,r

def subhalo_lens(n_draw):

	Sig_c = distance.critical_density(cosmopara,lenspara) # h M_sun/Mpc^2
	print Sig_c

	uni_m = np.random.rand(n_draw)
	#print uni_m
	m_d = MF.inverse_cdf(uni_m,ml,mu) # unit: M_sun

	Dl = distance.angular_distance(cosmopara,lenspara.zl)

	# kappa

	## b_sub (arc sec) & truncation radius

	b_com4 = (m_d/(np.pi*Sig_c))**2 # b_sub^3*b [Mpc^4]
	b_com4 = b_com4/Dl**4 # rad^4
	b_sub3 = b_com4/lenspara.b # b_sub^3 [rad^3]

	b_sub = (b_sub3)**(1.0/3.0) # rad
	b_sub = np.degrees(b_sub)*3600 # arc sec

	r_t = np.sqrt(b_sub*lenspara.b) # arc sec

	rt_m = np.radians(r_t/3600)*Dl  # Mpc
	Sig_d = m_d/(4.0*np.pi*rt_m**2)
	#print Sig_d
	k_d = Sig_d/Sig_c

	return m_d,k_d,b_sub,r_t
	

## -----------------

## ------ number of subhalo above a certain mass -------

# test subhalo total mass in subgroup [show in f_sub]

## draw 100

draw_s = 10000

x_s,y_s,rp = NFW_2d(draw_s) # centroid position [x,y] and distance to center

M_s,k_s,b_s,r_t = subhalo_lens(len(x_s))

'''
# positions within r_lim
x_s,y_s = x_s[rp<=2.0*lenspara.b],y_s[rp<=2.0*lenspara.b]

print len(x_s)

M_s,k_s,b_s,r_t = subhalo_lens(len(x_s))

print np.sum(2*k_s)
print k_s
#print b_s
print M_s
'''
plt.hist(np.log10(M_s),bins=100)
plt.show()

## Start to draw mass & r_t



