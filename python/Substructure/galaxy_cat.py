import numpy as np
import groupcat

basePath = '../../data/illustris_1'

ssNumber = 103 # snapshotnumber

mlow = 1e12
mhigh = 1e14

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])

list_end = list(GroupFirstSub).index(4294967295)

GroupFirstSub = GroupFirstSub[:list_end]


SubhaloMass = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloMass'])*1e10/cosmopara.h
GalaxyMass = SubhaloMass[GroupFirstSub]

SubhaloCM = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloCM'])
star_ms = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloMassType'])[:,4]*1e10/cosmopara.h

GalaxyCM = SubhaloCM[GroupFirstSub,:]
star_ms = star_ms[GroupFirstSub]

## use group mass as galaxy cretria 
mlow = 1e12
mhigh = 1e14

for i in range(SubhaloMass.size):
	if GalaxyMass[i]<mhigh:
		print i
		idx1 = i # idx_str
		break
for i in range(SubhaloMass.size):
	if GalaxyMass[i]<mlow:
		print i
		idx2 = i # idx_str
		break

GalaxyID = GroupFirstSub[idx1:idx2+1]
GalaxyMass = GalaxyMass[idx1:idx2+1]

#SubID_str = min(SubID)
#SubID_end = max(SubID)

#idx1 = list(GroupFirstSub).index(SubID_str)

GalaxyCM = GalaxyCM[idx1:idx2+1,:]
star_ms = star_ms[idx1:idx2+1]

catalog = open(basePath+'/Galaxy_'+str(ssNumber)+'.dat','w')

catalog.write('# [0]: Galaxy SubID w/ group mass between 10^'+str(np.log10(mlow))+'~10^'+str(np.log10(mhigh))+'\n')
catalog.write('# [1]-[3]: Galaxy CM in ckpc/h \n')
catalog.write('# [4]: Subhalo Mass in M_sun \n')
catalog.write('# [5]: Stellar mass in M_sun \n')

for i in range(GalaxyID.size):
	catalog.write(str(GalaxyID[i])+'    '+str(GalaxyCM[i,0])+'    '+str(GalaxyCM[i,1])+'    '+str(GalaxyCM[i,2])+'    '+str(GalaxyMass[i])+'    '+str(star_ms[i])+'\n')


