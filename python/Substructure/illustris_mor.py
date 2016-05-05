import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 103

Snapshot_num = 'snapshot_'+str(ssNumber)

## read in Galaxy catalog

catalog = basePath+'/Galaxy_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

# get morphology info [which band? u,g,i,H]

morFile = h5py.File(basePath+'/nonparametric_morphologies_z.hdf5')

# snapshot 103: z = 0.5

SubfindID = morFile['nonparmorphs']['snapshot_103']['SubfindID']
SubfindID = np.array(SubfindID)

## mask for the intersection
mask = [(i in GalaxyID) for i in SubfindID] ## This is important!! 
mask = np.array(mask)

Gini = morFile['nonparmorphs']['snapshot_103']['ACS-F606W']['CAMERA0']['GINI']
M20 = morFile['nonparmorphs']['snapshot_103']['ACS-F606W']['CAMERA0']['M20']

Gini = np.array(Gini)
M20 = np.array(M20)
#print Gini.size

Gini = Gini[mask]
print Gini[-10:]
M20 = M20[mask]
print M20[-10:]
SubfindID = SubfindID[mask]
#print Gini.size

def F(G,M20):

	#f = np.zeros(G.size)
	f = abs(-0.693*M20+4.95*G-3.85)
	mask = [G < (0.14*M20+0.778)]

	f[mask] = -1.0*f[mask]

	return f

F_index = F(Gini,M20)
print F_index.size


catalog = open(basePath+'/morphology_'+str(ssNumber)+'.dat','w')
catalog.write('# Galaxy SubID   Gini index   M20    F index'+'\n')

for i in range(SubfindID.size):
	catalog.write(str(SubfindID[i])+'    '+str(Gini[i])+'    '+str(M20[i])+'    '+str(F_index[i])+'\n')
'''

bulge = [F_index>0]
G_bulge = Gini[bulge]
M20_bulge = M20[bulge]

disk = [F_index<0]
G_disk = Gini[disk]
M20_disk = M20[disk]
print G_disk.size

plt.plot(M20_bulge,G_bulge,'+',label = 'bulge-dominated')
plt.plot(M20_disk,G_disk,'r+', label = 'disk-dominated')
plt.xlabel('M20')
plt.ylabel('Gini')
plt.xlim(0,-3)
plt.ylim(0.3,0.7)
plt.legend(loc = 4)
plt.show()


'''
