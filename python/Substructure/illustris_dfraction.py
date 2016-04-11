#import illustris_python as il
import numpy as np
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt

basePath = '../../data/illustris_3'

ssNumber = 99 # snapshotnumber
zl = 0.6
zs = 2.0

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27



GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])

mass_msun = np.zeros(len(GroupFirstSub))

#for i in range(len(GroupFirstSub)):

#	all_fields = groupcat.loadSingle(basePath,ssNumber,subhaloID = GroupFirstSub[i])
#	mass_msun[i] = all_fields['SubhaloMass']*1e10/cosmopara.h

SubhaloMass = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloMass'])

#print GroupFirstSub[0:10]
#print SubhaloMass[0:10]

j =0
for i in range(len(GroupFirstSub)):
	mass_msun[j] = SubhaloMass[i]*1e10/cosmopara.h
	j = j+1

mass_msun = mass_msun[mass_msun>1e9]

#plt.hist(mass_msun,bins = 10000)
#plt.xlim(1e8,1e12)
#plt.show()


''' 
#Einstein radius 
'''
theta_e = distance.EinsteinR(cosmopara,zl,zs,mass_msun)

plt.xlim(0,2.0)
plt.hist(theta_e,bins=100)
plt.show()
