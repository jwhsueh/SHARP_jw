""" This code test different profiles in 2D projected position """

import numpy as np
import NFWprofile as NFW
import scipy.stats

import matplotlib.pyplot as plt

## lenspara

class lenspara:
	q=0.3
	b=0.5 # arcsec

	zl=0.4
	zs=1.3

##------------

rs = NFW.scaleR(lenspara)
print rs
'''
""" NFW distirbution class """
class NFW_norm(scipy.stats.rv_continuous):
	def _pdf(self,r):
		x = r/rs
		return 1.0/x/(1+x)**2
		#return x

NFW_prob = NFW_norm(a=0,b=np.Inf,name = 'NFW_norm')

print NFW_prob.a, NFW_prob.b
check = NFW_prob.rvs(size=4)*10*rs

print check

#plt.hist(check)
#plt.show()
'''

# core radius
rc = 0.0 # arcsec

def NFW_profile(r):
	x = r/rs
	xc = rc/rs
	return 1.0/(x+xc)/(1+x)**2


x = np.random.rand(500000)*10*rs
check = NFW_profile(x)
plt.hist(check,bins = np.linspace(0,19,2000))
plt.xlim(0.1,1.0)
plt.ylim(0,2000)
plt.show()

'''

x = np.linspace(0,10*rs,10000)
y = NFW_profile(x)

plt.loglog(x,y)
plt.show()
'''