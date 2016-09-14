import numpy as np
import DistanceTool as distance
import matplotlib.pyplot as plt
import scipy.misc

class cosmopara:
	OM = 0.27
	h = 0.71

zl = np.array([0.2,0.6,1.0])

zs_0 = np.linspace(zl[0]+0.01,3.0,300)
zs_1 = np.linspace(zl[1]+0.01,3.0,300)
zs_2 = np.linspace(zl[2]+0.01,3.0,300)
#zs = 0.22

b=[]

def D(zl,zs):
	Dl = distance.angular_distance(cosmopara,zl)	
	Ds = distance.angular_distance(cosmopara,zs)
	Dls = Ds - Dl

	return Dls/Ds/Dl

diff=[]
D0,D1,D2 = [],[],[]

for i in range(len(zs_0)):
	D0.append(D(zl[0],zs_0[i]))
	D1.append(D(zl[1],zs_1[i]))
	D2.append(D(zl[2],zs_2[i]))

	#diff.append(scipy.misc.derivative(D,zs_se[i],dx=1e-6))

idx0 = D0.index(np.max(D0))
idx1 = D1.index(np.max(D1))
idx2 = D2.index(np.max(D2))

print zs_0[idx0],zs_1[idx1],zs_2[idx2]

plt.plot(zs_0,D0,'r',label='z=0.2')
plt.plot(zs_1,D0,'b',label='z=0.6')
plt.plot(zs_2,D0,'k',label='z=1.0')
plt.xlabel('source redshift')
plt.ylabel('$D_{LS}/D_L D_s$')
plt.legend(loc=4)
plt.show()