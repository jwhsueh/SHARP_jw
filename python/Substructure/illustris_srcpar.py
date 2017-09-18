import numpy as np
import matplotlib.pyplot as plt

path = '/Volumes/sting_1/snap99_'
#subID = 215763
subID = 179899
#re = 0.6864316# arcsec
re = 0.6373153
re = re*6.714 # kpc
a = 1.0/(1.0+0.6)
NN = 64

src_rad = 40 # pc
src_rad = 40/1000. # kpc

table = np.loadtxt(path+str(subID)+'/particle_'+str(subID)+'_dm.dat')
x,y,z = table[:,0]*a,table[:,1]*a,table[:,2]*a #Mpc

center_mass = np.array([np.sum(x)/len(x),np.sum(y)/len(x),np.sum(z)/len(x)])

x,y,z = x-center_mass[0],y-center_mass[1],z-center_mass[2]
x,y,z = x*1000.,y*1000.,z*1000. # kpc

## 0.5~1.5 re
r = np.sqrt(x**2+y**2) # proj1
#r = np.sqrt(x**2+z**2) # proj2
print np.min(r),np.max(r)
mask = np.logical_and(r>re-3*src_rad,r<re+3*src_rad) ## anulus w/ source rad
print mask[:10]
x,y,z = x[mask],y[mask],z[mask]

print len(x)
par_num = len(x)/(2.0*np.pi/6.0/src_rad*re)
print par_num
'''
smooth = np.zeros(len(x))
for i in range(len(x)):
	#print i
	xi,yi,zi = x[i],y[i],z[i]

	rd = np.sqrt((x-xi)**2+(y-yi)**2+(z-zi)**2)
	rd = np.sort(rd)
	sm = rd[NN+1]

	#print sm
	smooth[i] = sm

smooth = smooth[smooth<5.0]

plt.hist(smooth,bins=20)
plt.xlabel('smoothing length (kpc)')
plt.ylabel('number counts')
plt.title(str(subID))
#plt.show()
plt.savefig('smooth_pj_'+str(subID)+'.png')
'''