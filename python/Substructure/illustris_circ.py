import h5py
import numpy as np
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

basePath = '../../data/illustris_1'
ssNumber = 99

Snapshot_num = 'Snapshot_'+str(ssNumber)

circFile = h5py.File(basePath+'/stellar_circs.hdf5')

SubfindID_circ = circFile[Snapshot_num]['SubfindID']

CircAbove07Frac = circFile[Snapshot_num]['CircAbove07Frac']
beta = circFile[Snapshot_num]['CircTwiceBelow0Frac']
Js = circFile[Snapshot_num]['SpecificAngMom']

## read in Galaxy catalog

catalog = basePath+'/Galaxy_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
star_ms = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[4]) # stellar mass

GalaxyID_end = max(GalaxyID)

## mask for circFile
mask = [(i in GalaxyID) for i in SubfindID_circ] ## This is important!! 
mask = np.array(mask)

## star fraction& specific angular momentum

CircAbove07Frac = CircAbove07Frac[mask]
beta = beta[mask]
Js = Js[mask]
SubfindID_circ = SubfindID_circ[mask]

## indicator of specific angular momentum

Ms = []

for i in range(GalaxyID.size):
	if GalaxyID[i] in SubfindID_circ:
		Ms.append(star_ms[i])
	if i > GalaxyID_end: break

Ms = np.array(Ms)

AngMomIndex = np.log10(Js/Ms**(2./3.))

catalog = open(basePath+'/kinematics_'+str(ssNumber)+'.dat','w')
#catalog.write('# Galaxy SubID   Angular Momemtum   Disk star frac    Bulge star frac'+'\n')

catalog.write('# [0]: Galaxy SubID'+'\n')
catalog.write('# [1]-[2]: Disk star fraction & Bulge star fraction'+'\n')
catalog.write('# [3]: log10 SpecificAngMom/M_star^(2/3)'+'\n')


for i in range(SubfindID_circ.size):
	catalog.write(str(SubfindID_circ[i])+'    '+str(CircAbove07Frac[i])+'    '+str(beta[i])+'	'+str(AngMomIndex[i])+'\n')


'''

''' 
#Einstein radius 
'''
theta_e = distance.EinsteinR(cosmopara,zl,zs,sigma)
#print theta_e[0:10]


''' #select disk galaxy (different threshhold) 
'''

disk_e1 = theta_e[bStarFrac<0.6]
disk_ID1 = galaxyID[bStarFrac<0.6]

disk_e2 = theta_e[bStarFrac<0.4]
disk_ID2 = galaxyID[bStarFrac<0.4]

disk_e3 = theta_e[bStarFrac<0.2]
disk_ID3 = galaxyID[bStarFrac<0.2]

cat1 = open('beta_06.txt','w')
cat2 = open('beta_04.txt','w')
cat3 = open('beta_02.txt','w')

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
#print disk_hist1
#print disk_hist2
#print disk_hist3

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

plt.plot(np.linspace(0.,2.,40),fraction1,label = 'beta<0.6')
plt.plot(np.linspace(0.,2.,40),fraction2,'r',label = 'beta<0.4')
plt.plot(np.linspace(0.,2.,40),fraction3,'g--',label = 'beta<0.2')
plt.legend()
plt.xlim(0,1.2)
plt.xlabel('Einstein radius (arc sec)')
plt.ylabel('disk galaxy fraction')
plt.show()
'''
