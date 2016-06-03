import numpy as np
import snapshot
import h5py
import DistanceTool as distance

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = 99

class cosmopara:
	h = 0.704
	OM = 0.27

## read in Galaxy catalog

catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])

p_type = ['gas','dm','stars']

field = ['Masses','Velocities']

## read in snapshot header

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

## this function deal w/ galaxy on the boundary
def boundary(ci):
	
	#ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	for i in ci:
		if i<0.5*boxsize:
			i = i+boxsize

	return ci

## calculate angular momentum on each axis

angularM = np.zeros((GalaxyID.size,3))*np.nan

for i in range(GalaxyID.size):
