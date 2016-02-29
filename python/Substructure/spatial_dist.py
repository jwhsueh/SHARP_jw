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

##------------

rs = NFW.scaleR(lenspara)
r_end = 10*rs

''' CDF of aribitary distribution '''

ri = np.linspace(0,r_end,r_end*10000)
pdf_d = NFW.pdf(ri,rs) # discrete pdf
cdf_d = np.zeros(len(ri))

i = 1
while i < len(ri):
	cdf_d[i] = cdf_d[i-1]+pdf_d[i]*(ri[1]-ri[0])
	i = i+1

# normalization
cdf_d = cdf_d/max(cdf_d)

''' Inverse CDF '''

Ix = cdf_d
Iy = ri


rand = np.random.rand(10000)

check = np.interp(rand,Ix,Iy)

plt.hist(check, bins = np.linspace(0,r_end,100) )
plt.show()
