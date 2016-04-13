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


mlow = 1e12
mhigh = 1e14

GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])
GroupMass = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupMass'])*1e10/cosmopara.h

GroupFirstSub = GroupFirstSub[GroupMass>mlow]
GroupMass = GroupMass[GroupMass>mlow]

GroupFirstSub = GroupFirstSub[GroupMass<mhigh]

#SubhaloMass = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloMass'])
SubhaloVelDisp = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloVelDisp'])

#mass_msun = np.zeros(len(GroupFirstSub))
sigma = np.zeros(len(GroupFirstSub))

j = 0
for i in GroupFirstSub:

	#mass_msun[j] = SubhaloMass[i]*1e10/cosmopara.h
	sigma[j] = SubhaloVelDisp[i]
	#print sigma
	j = j+1





''' 
#Einstein radius 
'''
theta_e = distance.EinsteinR(cosmopara,zl,zs,sigma)
print theta_e[0:10]

plt.xlim(0,2.0)
plt.ylim(0,100)
plt.hist(theta_e)
plt.show()
