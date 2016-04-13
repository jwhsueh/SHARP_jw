import h5py
import numpy as np
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 99

mlow = 1e12
mhigh = 1e14

Snapshot_num = 'Snapshot_'+str(ssNumber)

circFile = h5py.File(basePath+'/stellar_circs.hdf5')

SubfindID = circFile[Snapshot_num]['SubfindID']

CircAbove07Frac = circFile[Snapshot_num]['CircAbove07Frac']




## cross check w/ GroupFirstSub

''' group catalog '''

zl = 0.6
zs = 2.0

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27


GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])
GroupMass = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupMass'])*1e10/cosmopara.h

galaxyID = GroupFirstSub[GroupMass>mlow]
GroupMass = GroupMass[GroupMass>mlow]

galaxyID = galaxyID[GroupMass<mhigh]

''' select of subhalo complete '''

dStarFrac = np.zeros(len(galaxyID))

j = 0

for i in range(len(SubfindID)):
	subID = SubfindID[i]
	#print subID
	if subID == galaxyID[j]:
		print subID
		dStarFrac[j] = CircAbove07Frac[i]
		j = j+1

		if j == (len(galaxyID)-1):
			break


SubhaloVelDisp = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloVelDisp'])


sigma = np.zeros(len(galaxyID))

j = 0
for i in galaxyID:
	sigma[j] = SubhaloVelDisp[i]
	j = j+1


''' 
#Einstein radius 
'''
theta_e = distance.EinsteinR(cosmopara,zl,zs,sigma)
#print theta_e[0:10]


''' select disk galaxy (f>0.4) '''

disk_e = theta_e[dStarFrac>0.4]
print len(disk_e)
print len(theta_e)

total_hist = np.histogram(theta_e,bins = np.linspace(0.,2.,41))
total_hist = total_hist[0]
#print total_hist

disk_hist = np.histogram(disk_e,bins = np.linspace(0.,2.,41))
disk_hist = disk_hist[0]
#print disk_hist

fraction = np.zeros(len(total_hist))
for i in range(len(total_hist)):
	if disk_hist[i]>0:
		#print disk_hist[i]
		fraction[i] = float(disk_hist[i])/float(total_hist[i])

#print fraction

plt.plot(np.linspace(0.,2.,40),fraction)
plt.show()
