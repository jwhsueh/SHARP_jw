import h5py
import numpy as np
import groupcat
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = 103

Snapshot_num = 'snapshot_'+str(ssNumber)

## read in Galaxy catalog

catalog = basePath+'/groups_'+str(ssNumber)+'/Galaxy_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog)

# get morphology info [which band? u,g,i,H]

morFile = h5py.File(basePath+'/nonparametric_morphologies_z.hdf5')

# snapshot 103: z = 0.5

SubfindID = morFile['nonparmorphs']['snapshot_103']['SubfindID']
SubfindID = np.array(SubfindID)
SubfindID = sorted(SubfindID)
print SubfindID[:50]

#mask = [(i in GalaxyID) for i in SubfindID]
#mask = np.array(mask)

#print SubfindID[mask]

#print np.array(SubfindID)
#data = morFile['nonparmorphs']['snapshot_103']['ACS-F606W']['CAMERA0']['GINI']
#print np.array(data)

Gini = morFile['nonparmorphs']['snapshot_103']['ACS-F606W']['CAMERA0']['GINI']

M20 = morFile['nonparmorphs']['snapshot_103']['ACS-F606W']['CAMERA0']['M20']

print set(GalaxyID).intersection(set(SubfindID))
'''
condition = []
for i in SubfindID:
	if i in GalaxyID:
		condition.append(True)
		print i
	else:
		condition.append(False)
'''
#print condition
#Gini = Gini[condition]
