import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

DanID = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'int',unpack=True, usecols=[0])
Dan_re = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'float',unpack=True, usecols=[1])

catalog = basePath+'/kinematics_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

bulge_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
disk_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'.dat'
Galaxy_ms = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[4])

Galaxy_V = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[6])
Galaxy_K = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[7])

## Re catalog
catalog_re = basePath+'/GalaxyRe_'+str(ssNumber)+'_x.dat'

Galaxy_re = np.loadtxt(catalog_re,dtype='float',unpack=True,usecols=[1])


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
Dan2_re = []

for i in range(GalaxyID.size):
	if GalaxyID[i] in DanID:
		Dan_ms.append(Galaxy_ms[i])
		bulge_dan.append(bulge_frac[i])
		disk_dan.append(disk_frac[i])
		Dan_V.append(Galaxy_V[i])
		Dan_K.append(Galaxy_K[i])

print len(Dan_ms)

for i in range(DanID.size):
	if DanID[i] in GalaxyID:
		Dan2_re.append(Dan_re[i])


## find kinematic info & mass for Dandan's selection


## mass

## disk criteria

cri = bulge_frac<0.6

DiskID = GalaxyID[cri]
print len(DiskID)

disk_cri = disk_frac[cri]
disk_ms = Galaxy_ms[cri]
disk_re = Galaxy_re[cri]
#disk_V = Galaxy_V[cri]
#disk_K = Galaxy_K[cri]

cri2 = disk_frac>0.4
disk_ms2 = Galaxy_ms[cri2]
disk_re2 = Galaxy_re[cri2]

bins = np.linspace(0,1.5,20)

plt.hist(Galaxy_re,bins,label = 'All galaxy',histtype='step')
plt.hist(disk_re,bins,color = 'g',label = 'bulge Star frac<0.6',histtype='step')
plt.hist(Dan2_re,bins,color = 'r', label='Dandan morphology pick',histtype='step')

plt.hist(disk_re2,bins,color = 'k',label = 'disk Star frac>0.4',histtype='step')
plt.legend()
plt.xlabel('Einstein radius (arc sec)')
plt.ylabel('galaxy count')
plt.title('Snapshot 99')
plt.show()

'''
plt.hist(disk_ms,bins,label = 'Bulge Star frac<0.6',)
plt.hist(Dan_ms,bins,color = 'r', label='Dandan morphology pick')
plt.hist(disk_ms2,bins,color = 'g',label = 'Disk Star frac>0.4',)
plt.legend()
plt.xlabel('M_sun')
plt.ylabel('galaxy count')
plt.title('Snapshot 99')
plt.show()
'''

'''

plt.ylim(30,16)
plt.plot(np.log10(Dan_ms),Dan_V,'+',label = 'morphology (Dandan)')
plt.plot(np.log10(disk_ms),disk_V,'r+',label = 'bulge star fraction<0.6')

#plt.plot(o_mass,o_frac,'o',mec ='g',mfc = 'none',label = 'overlap (beta<0.2)')
plt.xlabel('log10(Mass)')
plt.ylabel('V-band mag')
plt.legend()

plt.show()
'''