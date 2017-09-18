## This code create realization inside the lens halo w/ NFW profile

import numpy as np
import NFWprofile as NFW
import matplotlib.pyplot as plt
import sys

lens = 'B1422'
## critical density for lens
sig_c = 3.4e+10 #[M_Sun h^-1 arcsec^-2]
zl = 0.34

f_sub = float(sys.argv[1])
fold_name = sys.argv[2]
n_real = int(sys.argv[3])
n_start = int(sys.argv[4])

#path = './'+lens+'/'
path = '/Volumes/sting_1/subs/'
real_folder = 'real_'+fold_name+'/'

sublist_file = path+'sublist_'+lens+'_test.txt'

## subs fraction
# Here I want to draw from a uniform prior [0.001 to 0.05] first in log space
#fs_lo, fs_hi = np.log10(0.001), np.log10(0.05)
#f_sub = np.random.rand(n_real)*(fs_hi-fs_lo)+fs_lo
#f_sub = np.power(10,f_sub)
#print f_sub

## realization area
re = 0.75 # arcsec
#r = re*2
area = np.pi*(1.2*re-0.8*re)**2

## total mass for substructures
Mtot = f_sub*sig_c*area/2.0 # M_Sun h^-1
print Mtot

## --- read in sublist
#table = np.loadtxt('../../data/sub_gravlens/sublist_'+lens+'_test.txt')
table = np.loadtxt(sublist_file)
m_sub, prob = table[:,0], table[:,1]
r_sub, prob_r = table[:,2], table[:,3]

## --- save f_sub for each realization
real_f = open(path+real_folder+'real_fsub.txt','a')
real_n = open(path+real_folder+'real_nsub.txt','w')
#

for j in range(n_real):
	## draw mass
	## draw from probability
	n_inv = 5
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
		x_temp = np.random.rand(1)*3*re-1.5*re
		y_temp = np.random.rand(1)*3*re-1.5*re
		dis2_temp = np.sqrt(x_temp**2+y_temp**2)

		if np.logical_and(dis2_temp>(0.8*re),dis2_temp<1.2*re) :
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

	#plt.scatter(x_list,y_list,marker='o',facecolor='none',s=rt_list/np.min(rt_list))
	#plt.gca().set_aspect('equal')
	#plt.savefig('../../data/sub_gravlens/B1422_realization/real0.png')


#real0 = '../../data/sub_gravlens/B1422_realization/real0.txt'
#np.savetxt(real0,np.c_[m_list,x_list,y_list,bsub_list,rt_list])
