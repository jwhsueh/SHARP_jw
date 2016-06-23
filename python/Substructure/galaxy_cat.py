import numpy as np
import groupcat
import DistanceTool as distance

basePath = '../../data/illustris_1'

ssNumber = 85 # snapshotnumber
z = 1.0

## total mass
#mlow = 1e12
#mhigh = 1e14

#stellar mass
mlow = 5e9
mhigh = 10**11.8

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])

list_end = list(GroupFirstSub).index(4294967295)

GroupFirstSub = GroupFirstSub[:list_end]


SubhaloMass = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloMass'])*1e10/cosmopara.h
GalaxyMass = SubhaloMass[GroupFirstSub]

#GroupCM = groupcat.loadHalos(basePath,ssNumber, fields = ['GroupPos']) # Use GroupPos rather than CM!!! [next time to try subhaloPos]
GroupCM = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloPos'])
GroupCM = GroupCM[GroupFirstSub,:]
#GroupCM_x,GroupCM_y,GroupCM_z = GroupCM[:,0],GroupCM[:,1],GroupCM[:,2]

star_ms = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloMassType'])[:,4]*1e10/cosmopara.h
star_ms = star_ms[GroupFirstSub]

sigma = groupcat.loadSubhalos(basePath,ssNumber,fields = ['SubhaloVelDisp'])
sigma = sigma[GroupFirstSub]

sigma_h = 370.
sigma_l = 130.
'''
#photometry
V_mag = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloStellarPhotometrics'])[:,2]
K_mag = groupcat.loadSubhalos(basePath,ssNumber, fields = ['SubhaloStellarPhotometrics'])[:,3]

# apparent mag
DL = distance.luminosity_distance(cosmopara,z)*1e6 # pc
V_mag = V_mag+5.*(np.log10(DL)-1)
K_mag = K_mag+5.*(np.log10(DL)-1)
'''

#print star_ms

## rewrite this part!!! [str mass]
## use group mass as galaxy cretria 
## here!!
GalaxyID,GalaxyPos,Galaxy_ms,Galaxy_str,Galaxy_sig = [],[],[],[],[]

for i in range(sigma.size):
	if sigma[i]<sigma_h and sigma[i]>sigma_l:
		#print sigma[i]
		GalaxyID.append(GroupFirstSub[i])
		GalaxyPos.append(GroupCM[i,:])
		Galaxy_ms.append(GalaxyMass[i])
		Galaxy_str.append(star_ms[i])
		Galaxy_sig.append(sigma[i])


catalog = open(basePath+'/Galaxy_'+str(ssNumber)+'_sig.dat','w')

catalog.write('# [0]: Galaxy SubID w/ velocity dispertion between '+str(sigma_l)+'~'+str(sigma_h)+' km/s\n')
catalog.write('# [1]-[3]: Group pos in ckpc/h \n')
catalog.write('# [4]: Subhalo Mass in M_sun \n')
catalog.write('# [5]: Stellar mass in M_sun \n')
catalog.write('# [6]: 1d Velocity dispersion from all particles in km/s \n')
#catalog.write('# [6]-[7]: photometry V-band & K-band \n')

for i in range(len(GalaxyID)):
	catalog.write(str(GalaxyID[i])+'    '+str(GalaxyPos[i])[1:-1]+'    '+str(Galaxy_ms[i])+'    '+str(Galaxy_str[i])+'	'+str(Galaxy_sig[i])+'\n')
	#catalog.write(str(table[i,:])[1:-1]+'\n')

