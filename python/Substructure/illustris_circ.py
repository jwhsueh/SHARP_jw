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

gID = list(galaxyID)
sID = list(SubfindID)

sIndex = np.zeros(len(galaxyID))

j = 0
for element in gID:
	sIndex = sID.index(element)
	#print sIndex
	dStarFrac[j] = CircAbove07Frac[sIndex]
	j = j+1
	

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


''' select disk galaxy (different threshhold) '''

disk_e1 = theta_e[dStarFrac>0.2]
disk_ID1 = galaxyID[dStarFrac>0.2]

disk_e2 = theta_e[dStarFrac>0.4]
disk_ID2 = galaxyID[dStarFrac>0.4]

disk_e3 = theta_e[dStarFrac>0.1]
disk_ID3 = galaxyID[dStarFrac>0.1]

cat1 = open('CircAbove07Frac_02.txt','w')
cat2 = open('CircAbove07Frac_04.txt','w')
cat3 = open('CircAbove07Frac_01.txt','w')

files = [cat1,cat2,cat3]
IDs = [disk_ID1,disk_ID2,disk_ID3]
thetas = [disk_e1,disk_e2,disk_e3]

for i in range(len(files)):
	f = files[i]
	f.write('#SubfindID		theta_e(arc sec) \n')

	catID = IDs[i]
	cat_theta = thetas[i]

	for j in range(len(catID)):
		f.write(str(catID[j])+'		'+str(cat_theta[j])+'\n')

	f.close()


print len(disk_e1)
print len(disk_e2)
print len(disk_e3)
print len(theta_e)

total_hist = np.histogram(theta_e,bins = np.linspace(0.,2.,41))
total_hist = total_hist[0]
#print total_hist

disk_hist1 = np.histogram(disk_e1,bins = np.linspace(0.,2.,41))
disk_hist1 = disk_hist1[0]

disk_hist2 = np.histogram(disk_e2,bins = np.linspace(0.,2.,41))
disk_hist2 = disk_hist2[0]

disk_hist3 = np.histogram(disk_e3,bins = np.linspace(0.,2.,41))
disk_hist3 = disk_hist3[0]
#print disk_hist

fraction1 = np.zeros(len(total_hist))
fraction2 = np.zeros(len(total_hist))
fraction3 = np.zeros(len(total_hist))

for i in range(len(total_hist)):
	if disk_hist1[i]>0:
		#print disk_hist[i]
		fraction1[i] = float(disk_hist1[i])/float(total_hist[i])

	if disk_hist2[i]>0:
		#print disk_hist[i]
		fraction2[i] = float(disk_hist2[i])/float(total_hist[i])

	if disk_hist3[i]>0:
		#print disk_hist[i]
		fraction3[i] = float(disk_hist3[i])/float(total_hist[i])

#print fraction

plt.plot(np.linspace(0.,2.,40),fraction1,label = 'f>0.2')
plt.plot(np.linspace(0.,2.,40),fraction2,'r',label = 'f>0.4')
plt.plot(np.linspace(0.,2.,40),fraction3,'g--',label = 'f>0.1')
plt.legend()
plt.xlim(0,1.2)
plt.xlabel('Einstein radius (arc sec)')
plt.ylabel('disk fraction')
plt.show()

