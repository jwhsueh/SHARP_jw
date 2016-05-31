
import numpy as np
import DistanceTool as distance

""" Draw from a certain form of mass function (current: power-law) """

def pdf(m):
	a = -1.9 # power-law slope

	return m**a


""" Discrete CDF"""
def cdf_d(mi):

	pdf_d = pdf(mi) # discrete pdf
	cdf_d = np.zeros(len(mi))

	i = 1
	while i < len(mi):
		cdf_d[i] = cdf_d[i-1]+pdf_d[i]*(mi[1]-mi[0])
		i = i+1

	# normalization
	cdf_d = cdf_d/max(cdf_d)

	return cdf_d

""" Draw from inverse_cdf will get cloned distribution """

def inverse_cdf(ml,mu):

	mi = np.linspace(ml,mu,(mu-ml)/1e5)
	Ix = cdf_d(mi)
	Iy = mi

	#return np.interp(m,Ix,Iy)
	return Ix,Iy