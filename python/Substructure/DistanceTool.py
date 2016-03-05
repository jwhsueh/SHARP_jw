import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np

# read in cosmological parameters, redshifts (start & end)
# in Mpc
def angular_distance(cospara, z):

	h = cospara.h

	''' inverse E(z) '''
	#rEz= lambda x: (OM*(1.+x)**3+OK*(1+x)**2+OL)**(-0.5)
	rEz = lambda x: 1/Ez(cospara, x)


	# comoving distance 
	Dc = scipy.integrate.quad(rEz,0,z)
	Dc = 3000.0/h*Dc[0]

	DA = Dc/(1.0+z)

	return DA

def Ez(cospara, z):
	
	OM = cospara.OM
	OL = 1 - OM # flat universe
	OK = 1 - OM - OL

	return (OM*(1.+z)**3+OK*(1+z)**2+OL)**(0.5)


""" Critical Density Sigma_c in [h M_sun/Mpc^2] """

def critical_density(cospara, lenspara):

	c = 3e8 # m/s
	G = 6.67e-11 # m^3/kg/s^2

	## m to Mpc

	c = c/3.08e22  # Mpc/s
	G = G/(3.08e22)**3  # Mpc^3/kg/s^2

	zl,zs = lenspara.zl, lenspara.zs
	h = cospara.h

	Ds = angular_distance(cospara,zs)
	Dl = angular_distance(cospara,zl)
	Dls = Ds - Dl

	Sigma_c = c**2/(4.0*np.pi*G)*Ds/(Dl*Dls) # kg/Mpc^2
	Sigma_c = Sigma_c/2e30#/h # h M_sun/Mpc^2

	return Sigma_c

