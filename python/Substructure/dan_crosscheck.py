import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

DanID = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'int',unpack=True, usecols=[0])

catalog = basePath+'/kinematics_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

bulge_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
disk_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'.dat'
Galaxy_ms = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[4])

Galaxy_V = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[6])
Galaxy_K = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[7])

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

## Mass selection to Dandan's picl

Dan_ms = []
bulge_dan = []
disk_dan = []
Dan_V = []
Dan_K =[]
for i in range(GalaxyID.size):
	if GalaxyID[i] in DanID:
		Dan_ms.append(Galaxy_ms[i])
		bulge_dan.append(bulge_frac[i])
		disk_dan.append(disk_frac[i])
		Dan_V.append(Galaxy_V[i])
		Dan_K.append(Galaxy_K[i])
print len(Dan_ms)

## find kinematic info & mass for Dandan's selection


## mass

## disk criteria

cri = bulge_frac<0.6

DiskID = GalaxyID[cri]
print len(DiskID)

disk_cri = disk_frac[cri]
disk_ms = Galaxy_ms[cri]
disk_V = Galaxy_V[cri]
disk_K = Galaxy_K[cri]



plt.ylim(30,16)
plt.plot(np.log10(Dan_ms),Dan_V,'+',label = 'morphology (Dandan)')
plt.plot(np.log10(disk_ms),disk_V,'r+',label = 'bulge star fraction<0.6')

#plt.plot(o_mass,o_frac,'o',mec ='g',mfc = 'none',label = 'overlap (beta<0.2)')
plt.xlabel('log10(Mass)')
plt.ylabel('V-band mag')
plt.legend()

plt.show()