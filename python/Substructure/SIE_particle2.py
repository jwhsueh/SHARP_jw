import numpy as np
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt

## -----------

dm_mass = 6262734.72104
R = 1.8054702248995469/2
p_num = 884679
#e = 0.376
e = 0
q = 1-e

M_tot = p_num*dm_mass

out_file = open('/Volumes/sting_1/snap99_1798991/particle_1798991_dm.dat','w')

dm = np.loadtxt('/Volumes/sting_1/snap99_179899/particle_179899_dm.dat')

dm[:,0] = dm[:,0] - np.average(dm[:,0])
dm[:,1] = dm[:,1] - np.average(dm[:,1])
dm[:,2] = dm[:,2] - np.average(dm[:,2])

rx = np.append(dm[:,0],-1*dm[:,0])
rx = np.append(rx,-1*dm[:,0])
rx = np.append(rx,-1*dm[:,0])
rx = np.append(rx,dm[:,0])
rx = np.append(rx,dm[:,0])
rx = np.append(rx,dm[:,0])
rx = np.append(rx,-1*dm[:,0])

ry = np.append(dm[:,1],dm[:,1])
ry = np.append(ry,-1*dm[:,1])
ry = np.append(ry,-1*dm[:,1])
ry = np.append(ry,-1*dm[:,1])
ry = np.append(ry,-1*dm[:,1])
ry = np.append(ry,dm[:,1])
ry = np.append(ry,dm[:,1])

rz = np.append(dm[:,2],dm[:,2])
rz = np.append(rz,dm[:,2])
rz = np.append(rz,-1*dm[:,2])
rz = np.append(rz,dm[:,2])
rz = np.append(rz,-1*dm[:,2])
rz = np.append(rz,-1*dm[:,2])
rz = np.append(rz,-1*dm[:,2])

#H,xx,yy = np.histogram2d(rx,ry)
#plt.contour(np.log(H))
#plt.colorbar()
#plt.show()

print "start sampleing"

plist = 1./((1-e)*rx**2+(1+e)*ry**2+rz**2)
plist = plist/np.sum(plist) # normalized prob

p_idx = np.random.choice(plist.size,p_num, p=plist)
#print p_idx[:100]
print "sampeling done"

px,py,pz = [],[],[]

for idx in p_idx:
	px.append(rx[idx])
	py.append(ry[idx])
	pz.append(rz[idx])


px,py,pz = np.array(px),np.array(py),np.array(pz)

out_file.write('# nparticles '+str(px.size)+'\n')

for i in range(px.size):
	out_file.write(str(px[i])+'\t'+str(py[i])+'\t'+str(pz[i])+'\t'+str(dm_mass)+'\n')

out_file.close()

H,xx,yy = np.histogram2d(px,py)
plt.contour(np.log(H))
plt.colorbar()
plt.show()

#plt.scatter(par_x,par_y,s=1)
#plt.show()

#i = 0
#while (i <p_num):
#	cdf_i = np.random.random_sample()*(ru-rl)+rl

	# find r




