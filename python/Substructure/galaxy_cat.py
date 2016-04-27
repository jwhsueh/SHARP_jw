import numpy as np
import groupcat

basePath = '../../data/illustris_3'

ssNumber = 99 # snapshotnumber

mlow = 1e12
mhigh = 1e14

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])
GroupMass = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupMass'])*1e10/cosmopara.h

## use group mass as galaxy cretria 
mlow = 1e12
mhigh = 1e14

c1= GroupMass > mlow
c2 = GroupMass < mhigh

group1 = np.extract(c1, GroupFirstSub)
group2 = np.extract(c2, GroupFirstSub)

SubID = np.intersect1d(group1,group2)

catalog = open('Galaxy_099.dat','w')

catalog.write('# Galaxy SubID w/ group mass between 10^'+str(np.log10(mlow))+'~10^'+str(np.log10(mhigh))+'\n')

for i in SubID:
	catalog.write(str(i)+'\n')


