import numpy as np
import snapshot
import h5py
import DistanceTool as distance
import matplotlib.pyplot as plt

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = 99

class cosmopara:
	h = 0.704
	OM = 0.27

## read in Galaxy catalog
catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'_sig.dat'

GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])

p_type = ['gas','dm','stars']

field = ['Masses']

## read in snapshot header

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

#dm_ms =  header['MassTable'][1]*1e10/cosmopara.h # mass of dm particle [M_sun]
dm_ms =  header['MassTable'][1]*1e10

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

## projection function

def projection(c0,c1,c2,axis):

	if axis == 0: # y-z plane
		p0,p1 = c1,c2
	elif axis == 1: # x-z plane
		p0,p1 = c0,c2
	elif axis==2: # x-y plane
		p0,p1 = c0,c1

	return p0,p1

## -------start to plot particles

ID = 113727

idx = list(GalaxyID).index(ID)
print idx
GalaxyID = ID

subhalo_dm = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'dm')
coord = subhalo_dm['Coordinates'] # ckpc/h
dm_x,dm_y,dm_z = coord[:,0],coord[:,1],coord[:,2]

if (max(dm_x)-min(dm_x)> 0.5*boxsize): 
	print dm_x.size
	dm_x = boundary(dm_x)
	print 'cut x'
	print max(dm_x),min(dm_x)
	print dm_x.size
if (max(dm_y)-min(dm_y)> 0.5*boxsize): 
	print dm_y.size
	dm_y = boundary(dm_y)
	print 'cut y'
	print max(dm_y),min(dm_y)
	print dm_y.size
if (max(dm_z)-min(dm_z)> 0.5*boxsize): 
	print dm_z.size
	dm_z = boundary(dm_z)
	print 'cut z'
	print max(dm_z),min(dm_z)
	print dm_z.size

print 'dm coord loaded'

subhalo_gas = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'gas')
coord = subhalo_gas['Coordinates']
gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
	gas_x = boundary(gas_x)
if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
	gas_y = boundary(gas_y)
if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
	gas_z = boundary(gas_z)

print 'gas coord loaded'

subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'stars')
coord = subhalo_st['Coordinates']*a/1000. # Mpc/h
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

if (max(st_x)-min(st_x)> 0.5*boxsize): 
	gas_x = boundary(st_x)
if (max(st_y)-min(st_y)> 0.5*boxsize): 
	st_y = boundary(st_y)
if (max(st_z)-min(st_z)> 0.5*boxsize): 
	st_z = boundary(st_z)	

print 'stars coord loaded'

## change ref point to subhalo CM

dm_x,dm_y,dm_z = dm_x-CM_x[idx],dm_y-CM_y[idx],dm_z-CM_z[idx]
gas_x,gas_y,gas_z = gas_x-CM_x[idx],gas_y-CM_y[idx],gas_z-CM_z[idx]
st_x,st_y,st_z = st_x-CM_x[idx],st_y-CM_y[idx],st_z-CM_z[idx]
#print min(np.abs(dm_x))

'''
## projection & R_e
proj_axis = 2 #[2= in x ]

dm_p0,dm_p1 = projection(dm_x,dm_y,dm_z,proj_axis)
gas_p0,gas_p1 = projection(gas_x,gas_y,gas_z,proj_axis)
st_p0,st_p1 = projection(st_x,st_y,st_z,proj_axis)

print min(np.abs(dm_p0))
'''
# mass
#gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h
#st_ms = subhalo_st['Masses']*1e10/cosmopara.h

gas_ms = subhalo_gas['Masses']*1e10
st_ms = subhalo_st['Masses']*1e10


### write in mass file 

out_file = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(GalaxyID)+'.dat','w')

out_file.write('# nparticles '+str(dm_x.size+gas_x.size+st_x.size)+'\n')

for i in range(dm_x.size):
	out_file.write(str(dm_x[i])+' '+str(dm_y[i])+' '+str(dm_z[i])+' '+str(dm_ms)+'\n')

for i in range(gas_x.size):
	out_file.write(str(gas_x[i])+' '+str(gas_y[i])+' '+str(gas_z[i])+' '+str(gas_ms[i])+'\n')

for i in range(st_x.size):
	out_file.write(str(st_x[i])+' '+str(st_y[i])+' '+str(st_z[i])+' '+str(st_ms[i])+'\n')

out_file.close()
