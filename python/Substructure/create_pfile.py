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


p_type = ['gas','dm','stars']

field = ['Masses']

## read in snapshot header

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

dm_ms =  header['MassTable'][1]*1e10/cosmopara.h # mass of dm particle [M_sun]
#dm_ms =  header['MassTable'][1]*1e10

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor
boxsize = boxsize*a/cosmopara.h # kpc


## read in Galaxy catalog
catalog = '../../data/illustris_1/Galaxy_0'+str(snapNum)+'_test.dat'

GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])*a/1000./cosmopara.h # Mpc
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])*a/1000./cosmopara.h # Mpc
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])*a/1000./cosmopara.h # Mpc

#CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])*a/cosmopara.h # kpc
#CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])*a/cosmopara.h # kpc
#CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])*a/cosmopara.h # kpc


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

ID = 347491

idx = list(GalaxyID).index(ID)
print idx
GalaxyID = ID

subhalo_dm = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'dm')
coord = subhalo_dm['Coordinates']*a/1000.*cosmopara.h # Mpc
#coord = subhalo_dm['Coordinates']*a/cosmopara.h 
dm_x,dm_y,dm_z = coord[:,0],coord[:,1],coord[:,2]

if (max(dm_x)-min(dm_x)> 0.5*boxsize): 

	dm_x = boundary(dm_x)

if (max(dm_y)-min(dm_y)> 0.5*boxsize): 

	dm_y = boundary(dm_y)

if (max(dm_z)-min(dm_z)> 0.5*boxsize): 

	dm_z = boundary(dm_z)

print 'dm coord loaded'

subhalo_gas = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'gas')
coord = subhalo_gas['Coordinates']*a/1000.*cosmopara.h # Mpc
#coord = subhalo_gas['Coordinates']*a/cosmopara.h 
gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
	gas_x = boundary(gas_x)
if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
	gas_y = boundary(gas_y)
if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
	gas_z = boundary(gas_z)

print 'gas coord loaded'

subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'stars')
coord = subhalo_st['Coordinates']*a/1000.*cosmopara.h # Mpc
#coord = subhalo_st['Coordinates']*a/cosmopara.h 
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
## ----- mock for 5 kpc only ------- ##
r_cut = 0.01

dm_r = np.sqrt(dm_x**2+dm_y**2+dm_z**2)
gas_r = np.sqrt(gas_x**2+gas_y**2+gas_z**2)
st_r = np.sqrt(st_x**2+st_y**2+st_z**2)

dm_x,dm_y,dm_z = dm_x[dm_r<r_cut],dm_y[dm_r<r_cut],dm_z[dm_r<r_cut]
gas_x,gas_y,gas_z = gas_x[gas_r<r_cut],gas_y[gas_r<r_cut],gas_z[gas_r<r_cut]
st_x,st_y,st_z = st_x[st_r<r_cut],st_y[st_r<r_cut],st_z[st_r<r_cut]

# subset of star particles
#st_x,st_y,st_z = st_x[:70000],st_y[:70000],st_z[:70000]
'''
## change ref point to all positive
## ---- change ---- ##
ref_x,ref_y,ref_z = np.min(dm_x),np.min(dm_y),np.min(dm_z)
dm_x,dm_y,dm_z = dm_x-ref_x,dm_y-ref_y,dm_z-ref_z
gas_x,gas_y,gas_z = gas_x-ref_x,gas_y-ref_y,gas_z-ref_z
st_x,st_y,st_z = st_x-ref_x,st_y-ref_y,st_z-ref_z


# mass
gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h
st_ms = subhalo_st['Masses']*1e10/cosmopara.h


### write in mass file 

out_file = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(GalaxyID)+'_dm.dat','w')
out_file2 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(GalaxyID)+'_gas.dat','w')
out_file3 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(GalaxyID)+'_st.dat','w')

out_file.write('# nparticles '+str(dm_x.size)+'\n')
out_file2.write('# nparticles '+str(gas_x.size)+'\n')
out_file3.write('# nparticles '+str(st_x.size)+'\n')

for i in range(dm_x.size):
	out_file.write(str(dm_x[i])+'\t'+str(dm_y[i])+'\t'+str(dm_z[i])+'\t'+str(dm_ms)+'\n')

for i in range(gas_x.size):
	out_file2.write(str(gas_x[i])+'\t'+str(gas_y[i])+'\t'+str(gas_z[i])+'\t'+str(gas_ms[i])+'\n')

for i in range(st_x.size):
	out_file3.write(str(st_x[i])+'\t'+str(st_y[i])+'\t'+str(st_z[i])+'\t'+str(st_ms[i])+'\n')

out_file.close()
