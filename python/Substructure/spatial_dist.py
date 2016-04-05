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

#rs = NFW.scaleR(cosmopara,lenspara)
'''use set_halopara to replace'''
halopara = NFW.set_halopara(cosmopara,lenspara)
'''
print halopara.M_200
print halopara.c_200
print halopara.rs
print halopara.r_200
print halopara.rho_s
'''
rs = halopara.rs[0] # arc sec
print rs
r_end = 100*lenspara.b

## -----substructure setting---------##

ml,mu = 1e6, 1e8 # lower & upper limits
f_sub = 0.01 # substructure mass fraction

## -------------

def NFW_2d(n_draw):
#n_draw = 100000

	uni_r = np.random.rand(n_draw)

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

	return x_d,y_d,rp_d,r_d


def uniform_2d(n_draw):
	x = lenspara.xc+np.random.uniform(x_lim[0],x_lim[1],n_draw)
	y = lenspara.yc+np.random.uniform(y_lim[0],y_lim[1],n_draw)

	r = np.sqrt(x_d**2+y_d**2)

	return x,y,r

def subhalo_lens(r_d):

	n_draw = len(r_d)

	Sig_c = distance.critical_density(cosmopara,lenspara) # h M_sun/Mpc^2
	print Sig_c

	uni_m = np.random.rand(n_draw)
	#print uni_m
	m_d = MF.inverse_cdf(uni_m,ml,mu) # unit: M_sun

	Dl = distance.angular_distance(cosmopara,lenspara.zl)

	# kappa

	## b_sub (arc sec) & truncation radius

	r_t = r_d*(m_d/NFW.enclose_mass(cosmopara,r_d,lenspara.zl,halopara))**(1.0/3.0) # arc sec

	rt_M = distance.arcs2mpc(cosmopara,r_t,lenspara.zl)  # Mpc
	Sig_d = m_d/(np.pi*rt_M**2)/cosmopara.h # h M_sun/Mpc^2
	#print Sig_d

	b_sub = m_d/(np.pi*rt_M*Sig_c*cosmopara.h) # Mpc
	b_sub = distance.mpc2arcs(cosmopara,b_sub,lenspara.zl) # arc sec

	k_d = Sig_d/Sig_c

	return m_d,k_d,b_sub,r_t
	

## -----------------

## ------ number of subhalo above a certain mass -------

# test subhalo total mass in subgroup [show in f_sub]

## draw 100

draw_s = 100

x_s,y_s,rp,rd = NFW_2d(draw_s) # centroid position [x,y] projected distance to center [rp], and actual distance [rd]

#M_s,k_s,b_s,r_t = subhalo_lens(len(x_s))


# positions within r_lim
x_s,y_s,rd = x_s[rp<=2.0*lenspara.b],y_s[rp<=2.0*lenspara.b],rd[rp<=2.0*lenspara.b]

print len(rd)

M_s,k_s,b_s,r_t = subhalo_lens(rd)

print np.sum(2*k_s)

#print k_s
#print b_s
#print r_t
#print M_s


plt.scatter(x_s,y_s)

th = np.linspace(0,2*np.pi,100)
#plt.plot(0.5*np.cos(th),0.5*np.sin(th),'r', linewidth = 1.0)
plt.plot(0.55*np.cos(th),0.55*np.sin(th),'r-', linewidth = 1.0)
plt.plot(0.45*np.cos(th),0.45*np.sin(th),'r-', linewidth = 1.0)

#plt.hist(np.log10(M_s),bins=100)
#plt.hist(rp,bins = 50)
#plt.xlim(0,1.)

plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('NFW 2D realization')
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.axis('equal')
#plt.savefig('NFW_2D_realization.png')

plt.show()

## Start to draw mass & r_t



