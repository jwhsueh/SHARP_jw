import numpy as np
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt

## -----------

dm_mass = 6262734.72104/2
R = 1.8054702248995469/2
p_num = 884679
e = 0.235
q = 1-e

M_tot = p_num*dm_mass

out_file = open('/Volumes/sting_1/snap99_1798992/particle_1798992_dm.dat','w')

## -----------
# pdf(r) is r^-1
idv = int(1e2+1)
rb = np.linspace(7e-2,R,idv)
print rb
rb2 = np.linspace(0,7e-2,int(idv/(R)*7e-2))
A_tot = np.log(R)-np.log(rb[0])
#A_tot = 1./rb[0]-1./rb[-1]

b_num = np.zeros(idv-1)

for i in range(idv-1):
	area = np.log(rb[i+1])-np.log(rb[i])
	#area = 1./rb[i]-1./rb[i+1]
	frac = area/A_tot
	#print frac
	b_num[i] = np.round(p_num*frac)

b_num = b_num.astype(int)/2
#print np.max(b_num)

b_num_core = np.zeros(rb2.size-1)
b_num_core.fill(np.max(b_num))
#print rb2[1:-1]

b_num = np.append(b_num_core[:-1],b_num)
rb = np.append(rb2[1:-1],rb)
#print b_num.size,rb.size

#plt.plot(rb[:-1],b_num)
#plt.plot(rb2[1:],b_num_core)
#plt.show()

## re-normalize

b_num = np.round(b_num*p_num/np.sum(b_num)).astype(int)
#print np.sum(b_num)

## ---- randomly generate coordinate

par_x,par_y,par_z = [],[],[]

for i in range(rb.size-1):
	pn = 0
	idr = rb[i+1]-rb[i]
	rint,rend = rb[i],rb[i+1]

	print "#### bin "+str(i)
	print "## b_num = "+str(b_num[i])
	print "####"
	while(pn<=b_num[i]):
		# draw r randomly w/i the bin
		ri = np.random.random_sample()*idr+rint
		# draw angles
		theta = np.random.random_sample()*2*np.pi
		phi = np.random.random_sample()*np.pi

		# coordinate
		par_x.append(ri*np.cos(theta)*np.sin(phi))
		par_y.append(ri*np.sin(theta)*np.sin(phi))
		par_z.append(ri*np.cos(phi))

		pn=pn+1

print par_x[:10], len(par_x)
hp_num = len(par_x)

par_x,par_y,par_z = np.array(par_x),np.array(par_y),np.array(par_z)

# bring in ellipticity
par_x = par_x*(1-e)

out_file.write('# nparticles '+str(hp_num)+'\n')

for i in range(hp_num):
	out_file.write(str(par_x[i])+'\t'+str(par_y[i])+'\t'+str(par_z[i])+'\t'+str(dm_mass)+'\n')

out_file.close()

#plt.scatter(par_x,par_y,s=1)
#plt.show()

#i = 0
#while (i <p_num):
#	cdf_i = np.random.random_sample()*(ru-rl)+rl

	# find r




