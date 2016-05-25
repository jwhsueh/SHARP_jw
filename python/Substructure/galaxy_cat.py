import numpy as np
import groupcat
import DistanceTool as distance

basePath = '../../data/illustris_1'

ssNumber = 108 # snapshotnumber
z = 0.4

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

GroupCM = groupcat.loadHalos(basePath,ssNumber, fields = ['GroupPos'])

star_ms = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloMassType'])[:,4]*1e10/cosmopara.h
star_ms = star_ms[GroupFirstSub]

#photometry
V_mag = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloStellarPhotometrics'])[:,2]
K_mag = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloStellarPhotometrics'])[:,3]

# apparent mag
DL = distance.luminosity_distance(cosmopara,z)*1e6 # pc
V_mag = V_mag+5.*(np.log10(DL)-1)
K_mag = K_mag+5.*(np.log10(DL)-1)

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
		idx2 = i # idx_end
		break

GalaxyID = GroupFirstSub[idx1:idx2+1]
GalaxyMass = GalaxyMass[idx1:idx2+1]

#SubID_str = min(SubID)
#SubID_end = max(SubID)

#idx1 = list(GroupFirstSub).index(SubID_str)

GalaxyCM = GroupCM[idx1:idx2+1,:]
star_ms = star_ms[idx1:idx2+1]

catalog = open(basePath+'/Galaxy_'+str(ssNumber)+'.dat','w')

catalog.write('# [0]: Galaxy SubID w/ group mass between 10^'+str(np.log10(mlow))+'~10^'+str(np.log10(mhigh))+'\n')
catalog.write('# [1]-[3]: Group pos in ckpc/h \n')
catalog.write('# [4]: Subhalo Mass in M_sun \n')
catalog.write('# [5]: Stellar mass in M_sun \n')
catalog.write('# [6]-[7]: photometry V-band & K-band \n')

for i in range(GalaxyID.size):
	catalog.write(str(GalaxyID[i])+'    '+str(GalaxyCM[i,0])+'    '+str(GalaxyCM[i,1])+'    '+str(GalaxyCM[i,2])+'    '+str(GalaxyMass[i])+'    '+str(star_ms[i])+'    '+str(V_mag[i])+'    '+str(K_mag[i])+'\n')


