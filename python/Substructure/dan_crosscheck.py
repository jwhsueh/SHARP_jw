import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 99

mlow = 1e12
mhigh = 1e14

Dan = np.loadtxt('JenWei_table_disk_in_x.dat')
Frac = np.loadtxt('CircAbove07Frac_04.txt')

mID = Dan[:,0] # morphology ID
fID = Frac[:,0] # dynamics ID

mID = list(mID)
fID = list(fID)

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])
GroupMass = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupMass'])*1e10/cosmopara.h

#galaxyID = GroupFirstSub[GroupMass>mlow]
#GroupMass = GroupMass[GroupMass>mlow]

galaxyID = GroupFirstSub

''' select of subhalo complete '''

''' Disk star fraction '''

Snapshot_num = 'Snapshot_'+str(ssNumber)

circFile = h5py.File(basePath+'/stellar_circs.hdf5')

SubfindID = circFile[Snapshot_num]['SubfindID']

CircAbove07Frac = circFile[Snapshot_num]['CircAbove07Frac']

m_frac = np.zeros(len(mID))
f_frac = np.zeros(len(fID))

sID = list(SubfindID)

''' Subhalo mass '''

SubhaloMass = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloMass'])*1e10/cosmopara.h

gID = list(galaxyID)

m_mass = np.zeros(len(mID))
f_mass = np.zeros(len(fID))

j = 0
for element in mID:
	mIndex = gID.index(element)
	mIndex2 = sID.index(element)
	#print sIndex
	m_mass[j] = SubhaloMass[mIndex]
	m_frac[j] = CircAbove07Frac[mIndex2]
	j = j+1

j = 0
for element in fID:
	fIndex = gID.index(element)
	fIndex2 = sID.index(element)
	#print sIndex
	f_mass[j] = SubhaloMass[fIndex]
	f_frac[j] = CircAbove07Frac[fIndex2]
	j = j+1	

#m_mass = m_mass[m_mass>1e12]
#m_mass = m_mass[m_mass<1e14]
#m_frac = m_frac[m_mass>1e12]
#m_frac = m_frac[m_mass<1e14]

m_mass = np.array(np.log10(m_mass))
f_mass = np.log10(f_mass)

print m_mass[0:10]
print m_frac[0:10]

#plt.xlim(12.,14.)
plt.scatter(m_mass[0:10],m_frac[0:10],'r')#,label = 'morphology')
#plt.scatter(f_mass,f_frac)#,label = 'f>0.4')
plt.xlabel('log10(Mass)')
plt.ylabel('Disk star fraction')
#plt.legend()

plt.show()