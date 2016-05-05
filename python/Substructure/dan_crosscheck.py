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

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

## Mass selection to Dandan's picl

Dan_ms = []
bulge_dan = []
disk_dan = []
for i in range(GalaxyID.size):
	if GalaxyID[i] in DanID:
		Dan_ms.append(Galaxy_ms[i])
		bulge_dan.append(bulge_frac[i])
		disk_dan.append(disk_frac[i])
print len(Dan_ms)

## find kinematic info & mass for Dandan's selection


## mass

## disk criteria

cri = bulge_frac<0.6

DiskID = GalaxyID[cri]
print len(DiskID)

disk_cri = bulge_frac[cri]
disk_ms = Galaxy_ms[cri]


#plt.xlim(1e12,1e14)
plt.plot(Dan_ms,bulge_dan,'+',label = 'morphology pick')
plt.plot(disk_ms,disk_cri,'r+',label = 'bulge star fraction<0.6')
#plt.plot(o_mass,o_frac,'o',mec ='g',mfc = 'none',label = 'overlap (beta<0.2)')
plt.xlabel('log10(Mass)')
plt.ylabel('Disk star fraction (beta)')
plt.legend()

plt.show()