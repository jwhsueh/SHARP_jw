import numpy as np

## ------- This code create a list of ramdom mass and a list of distance to halo centre
## ------- from a certain mass range and distance range
## ------- calculate corresponding probability from the mass function [dn/dm ~ M^-1.9]
## ------- radial distribution [dn/dr ~ M^-1]

## mass limit
m0,m1 = 6,9 # log10(M)

## mass function
def p_mass(mass):
	prob = mass**(-1.9)

	return prob
'''
## radial limit
re = 0.75 # arcsec
r0,r1 = 2*re,200*re
r0,r1 = np.log10(r0), np.log10(r1)

## radial function
def p_radial(distance):
	prob = np.exp(-1.0*distance**0.678) ## Einasto profile

	return prob
'''
## number of mass
n = 1e6

## draw from log10 space
log_m = np.random.rand(n)*(m1-m0)+m0
mass = np.power(10,log_m)
mass = mass[mass>10**m0]

## probability 
prob = p_mass(mass)
prob = prob/np.sum(prob) # normalize

'''
## draw distance from log10 space
log_r = np.random.rand(n)*(r1-r0)+r0
distance = np.power(10,log_r)

prob_r = p_radial(distance)
prob_r = prob_r/np.sum(prob_r)
'''
#print prob[:10]
#print mass[:10]

sublist = '../../data/sub_gravlens/sublist_CDM.txt'
np.savetxt(sublist,np.c_[mass,prob])