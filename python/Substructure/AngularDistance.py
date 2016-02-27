import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np

# read in cosmological parameters, redshifts (start & end)
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