import numpy as np
import pjaffe as pj
import matplotlib.pyplot as plt
import queue, time, threading, datetime, sys

lens = 'B1422'
## critical density for lens
sig_c = 3.4e+10 #[M_Sun h^-1 arcsec^-2]

fsub = float(sys.arg[1])
n_real = int(sys.arg[2])
n_start = int(sys.arg[3])

sublist_file = './'+lens+'/sublist_'+lens+'.txt'

## subs fraction
# Here I want to draw from a uniform prior [0.001 to 0.05] first in log space
#fs_lo, fs_hi = np.log10(0.001), np.log10(0.05)
#f_sub = np.random.rand(n_real)*(fs_hi-fs_lo)+fs_lo
#f_sub = np.power(10,f_sub)
#print f_sub

## realization area
re = 0.75 # arcsec
r = re*2
area = np.pi*r**2

## total mass for substructures
Mtot = f_sub*sig_c*area/2.0 # M_Sun h^-1
#print Mtot

## --- read in sublist
#table = np.loadtxt('../../data/sub_gravlens/sublist_'+lens+'_test.txt')
table = np.loadtxt(sublist_file)
m_sub, prob = table[:,0], table[:,1]
r_sub, prob_r = table[:,2], table[:,3]

## --- save f_sub for each realization
real_f = open('../../data/sub_gravlens/'+lens+'_realization1/real_fsub.txt','a')
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
		
		if np.sum(m_list)>Mtot[j]:
			break

	print len(m_list)
	f_dif = (np.sum(m_list) - Mtot[j])/(sig_c*area)*2.0
	print f_dif
	f_sub_j = np.sum(m_list)/(sig_c*area)*2.0
	real_f.write(str(f_sub_j)+'\n')
	#f_sub[j] = f_sub_j
	#print r_list[:10]

	## draw spatial distribution
	## ----- 2D distribution is uniform
	## ----- draw 2D position first and then draw a distance rsub > r_2d

	x_list,y_list,r_list = np.empty(0),np.empty(0),np.empty(0)
	i=0
	while i<len(m_list):
		#print "####"
		#print i
		#print '####'
		x_temp = np.random.rand(1)*2*r-r
		y_temp = np.random.rand(1)*2*r-r
		dis2_temp = x_temp**2+y_temp**2

		if (dis2_temp) < r**2:
			#print dis2_temp, x_temp,y_temp
			x_list = np.append(x_list,x_temp)
			y_list = np.append(y_list,y_temp)
			i=i+1
		## draw rsub
			r_temp=0.0
			while r_temp<dis2_temp:
				r_temp = np.random.choice(r_sub,1,p=prob_r)

			r_list = np.append(r_list,r_temp)



	## calculate bsub, r_t
	bsub_list = pj.bsub(m_list,r_list,sig_c,re)
	rt_list = pj.rt(r_list,bsub_list,re)

	real_i = '../../data/sub_gravlens/'+lens+'_realization1/real'+str(j+n_start)+'.txt'
	np.savetxt(real_i,np.c_[m_list,x_list,y_list,bsub_list,rt_list])

	#plt.scatter(x_list,y_list,marker='o',facecolor='none',s=rt_list/np.min(rt_list))
	#plt.gca().set_aspect('equal')
	#plt.savefig('../../data/sub_gravlens/B1422_realization/real0.png')


#real0 = '../../data/sub_gravlens/B1422_realization/real0.txt'
#np.savetxt(real0,np.c_[m_list,x_list,y_list,bsub_list,rt_list])
