import numpy as np

def bsub(Msub,rs,sig_c,b):
	return (Msub/sig_c)**(2./3.)*(3.0*b)**(1./3.)/np.pi/rs

def rt(rs,bsub,b):
	return rs**(3./2.)/np.sqrt(b)*np.sqrt(np.pi*bsub/3.0/b)