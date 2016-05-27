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
ssNumber = 108

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

