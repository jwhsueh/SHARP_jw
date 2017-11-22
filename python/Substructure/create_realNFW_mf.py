## This code create realization inside the lens halo w/ NFW profile

import numpy as np
import NFWprofile as NFW
import sys

lens = 'MG0414'
## critical density for lens
sig_c = 3.4e+10 #[M_Sun h^-1 arcsec^-2] here!!!
zl = 0.96

f_sub = float(sys.argv[1])
fold_name = sys.argv[2]
n_real = int(sys.argv[3])
n_start = int(sys.argv[4])

#path = './'+lens+'/'
path = '/Volumes/sting_1/subs/'+lens
real_folder = '/real_'+fold_name+'/'

## realization area
re = 1.11 # arcsec
#r = re*2
f_low,f_hi = 0.8,1.5
area = np.pi*(f_hi*re-f_low*re)**2

## total mass for substructures
Mtot = f_sub*sig_c*area/2.0 # M_Sun h^-1
print Mtot

m_sub = np.linspace(6,9,1000)
m_sub = 10**m_sub

m_c = 0.0

def massf(mass,m_c):
	dn = mass**(-1.9)*(1+m_c/mass)**(-1.3)

	return dn

## --- save f_sub for each realization
real_f = open(path+real_folder+'real_fsub.txt','w')
real_n = open(path+real_folder+'real_nsub.txt','w')

## probability of substurcture mass
prob = massf(m_sub,m_c)
prob = prob/np.sum(prob)

for j in range(n_real):
	## draw mass
	## draw from probability
	n_inv = 1
	m_list = np.empty(0)
	#r_list = np.empty(0)

	print '####'
	print 'realization '+str(j+n_start)
	print '####'

	#n = 0
	#m_sum = 0
	while True:
		inv_list = np.random.choice(m_sub,n_inv,p=prob)
		#inv_rlist = np.random.choice(r_sub,n_inv,p=prob_r)
		m_list = np.append(m_list,inv_list)
		#r_list = np.append(r_list,inv_rlist)
		#print np.sum(m_list)
		
		if np.sum(m_list)>Mtot:
			break

	print len(m_list)
	real_n.write(str(len(m_list))+'\n')

	f_dif = (np.sum(m_list) - Mtot)/(sig_c*area)*2.0
	print f_dif
	f_sub_j = np.sum(m_list)/(sig_c*area)*2.0
	real_f.write(str(f_sub_j)+'\n')
	#f_sub[j] = f_sub_j
	#print r_list[:10]

	## draw spatial distribution
	## ----- 2D distribution is uniform

	x_list,y_list = np.empty(0),np.empty(0)
	i=0
	while i<len(m_list):
		#print "####"
		#print i
		#print '####'
		x_temp = np.random.rand(1)*2.0*f_hi*re-f_hi*re
		y_temp = np.random.rand(1)*2.0*f_hi*re-f_hi*re
		dis2_temp = np.sqrt(x_temp**2+y_temp**2)

		if np.logical_and(dis2_temp>(f_low*re),dis2_temp<f_hi*re) :
			#print dis2_temp, x_temp,y_temp
			x_list = np.append(x_list,x_temp)
			y_list = np.append(y_list,y_temp)
			i=i+1



	## calculate ks, rs
	r200_list = NFW.r200(zl,m_list)
	c200_list = NFW.c200(zl,m_list)
	rs_list = r200_list/c200_list

	#rhos_list = NFW.rho_s(m_list,r200_list,c200_list)
	#ks_list = rhos_list*rs_list/sig_c


	real_i = path+real_folder+'real'+str(j+n_start)+'.txt'
	np.savetxt(real_i,np.c_[m_list,x_list,y_list,r200_list,c200_list])

