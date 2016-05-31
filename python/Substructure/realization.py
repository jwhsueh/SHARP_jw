""" This code generate realizations under fixed f_sub and call gravlens findimg """

import numpy as np
import matplotlib.pyplot as plt

import MassFunc as MF
import DistanceTool as distance
import Basic_class as bClass

## initiate lens
lens_name = 'B1422'
lensPath = '../../models/lens_info/'

# Lens setup file
lens_setup = np.loadtxt(lensPath+lens_name+'_setup.dat')
lens_obs = np.loadtxt(lensPath+lens_name+'_obs.dat')

#B1422 = bClass.Lens(z_lens= 0.536,z_src = 3.62,ellipticity = 0.4, Re = 0.7511, x_cen = 0.74,y_cen = -0.66)
Lens = bClass.Lens(lens_setup)
Lens.obsData(lens_obs)

## realization boundary
x_lim = np.array([Lens.xc-2.*Lens.b,Lens.xc+2.*Lens.b])
y_lim = np.array([Lens.yc-2.*Lens.b,Lens.yc+2.*Lens.b])


## substructure mass fraction
f_sub = 0.01

sigma_c = Lens.critical_density()*cospara.h # M_sun/Mpc^2
kappa_tot = f_sub*2.0/sigma_c # M_sun/Mpc^2

## cosmology parameter
cospara = bClass.Cosmology()

## 2D position + mass
## truncation radius (fixed at 1/60 -- F&K)

def uniform_2d(n_draw):
	x = np.random.uniform(x_lim[0],x_lim[1],n_draw)
	y = np.random.uniform(y_lim[0],y_lim[1],n_draw)

	dx = x - Lens.xc
	dy = y - Lens.yc

	r = np.sqrt(dx**2+dy**2)

	x = x[r<=2.*Lens.b]
	y = y[r<=2.*Lens.b]
	#r = r[r<=2.*Lens.b]

	return x,y

def get_mass(n_draw):
	draws = np.random.rand(n_draw)

	m_min,m_max = 1e6, 1e8 # lower & upper limits

	MF_x,MF_y = MF.inverse_cdf(m_min,m_max)

	m = np.interp(draws,MF_x,MF_y)

	return m

def m2b_jaffe(m_tot):
	# M_sub = pi*b*b_sub^3/2 * sigma_c

	b_mpc = distance.arcs2mpc(Lens.b)

	bs_15 = m_tot/np.pi/b_mpc/sigma_c
	bs_mpc = np.power(bs_15,2./3.)
	bs = distance.mpc2arcs(bs_mpc)

	return bs

def jaffe_ks(b_sub,rt):

	K = lambda r: b_sub/2.0*(1./r-1./np.sqrt(r**2+rt**2))

	m_in = 2.0*np.pi*np.quad(0,b_sub,K)[0]

	ks = m_in/np.pi/b_sub**2

	return ks


## Function: create realization

def set_realization(): # need to add realization number

##---- Here we start the drawing ----##

	unit_draw = 10 # 1 unit for the drawing

	# position
	xi,yi = uniform_2d(unit_draw)

	unit_sub = xi.size # number of substructure in this unit draw

	mi = get_mass(unit_sub)

	# mass to pseudo-Jaffe profile
	b_i = m2b_jaffe(mi)
	rt_i = np.sqrt(b_i*Lens.b)

	# pjaffe ks calculation
	ks_sub = jaffe_ks(b_i,rt_i)

	## set up counting approach


