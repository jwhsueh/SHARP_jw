import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

DanID = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'int',unpack=True, usecols=[0])
DanID.sort()
#Dan_re = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'float',unpack=True, usecols=[1])

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

## Dandan's full catalog

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'.dat'
LensID = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Re = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])
Mass = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[1])/0.704


''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

## Mass selection to Dandan's pick
## also find DMfrac

## now ID is sorted

id_i = 11625
id_e = 207298

DanID_mc = np.intersect1d(DanID[DanID>=id_i],DanID[DanID<=id_e]) # mass cut

# Re
DanRe_mc = []
DanMass_mc = []
DanDMfrac_mc = []
for i in range(LensID.size):
	if LensID[i] in DanID_mc:
		DanRe_mc.append(Re[i])
		DanMass_mc.append(Mass[i])
		DanDMfrac_mc.append(DMfrac[i])

## Galaxy Re & DM frac

# use index

Galaxy_re = np.zeros(GalaxyID.size)*np.nan
Galaxy_DMfrac = np.zeros(GalaxyID.size)*np.nan

for i in range(GalaxyID.size):
	if GalaxyID[i] in LensID:
		idx = list(LensID).index(GalaxyID[i])
		Galaxy_re[i] = Re[idx]
		Galaxy_DMfrac[i] = DMfrac[idx]

print Galaxy_DMfrac.size, Galaxy_re.size
#####

cri = bulge_frac<0.6

DiskID = GalaxyID[cri]
print len(DiskID)

disk_cri = disk_frac[cri]
disk_ms = Galaxy_ms[cri]
disk_re = Galaxy_re[cri]
disk_DMfrac = Galaxy_DMfrac[cri]
#disk_V = Galaxy_V[cri]
#disk_K = Galaxy_K[cri]

cri2 = disk_frac>0.4
disk_ms2 = Galaxy_ms[cri2]
disk_re2 = Galaxy_re[cri2]
disk_DMfrac2 = Galaxy_DMfrac[cri2]

## find re and mass for DiskID

bins = np.linspace(0.2,1.0,20)

plt.hist(Galaxy_DMfrac,bins,label = 'All galaxy',histtype='step')
plt.hist(disk_DMfrac,bins,color = 'g',label = 'bulge Star frac<0.6',histtype='step')
plt.hist(DanDMfrac_mc,bins,color = 'r', label='Dandan morphology pick',histtype='step')

plt.hist(disk_DMfrac2,bins,color = 'k',label = 'disk Star frac>0.4',histtype='step')
plt.legend(loc =1)
plt.xlabel('DM frac w/i Re')
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