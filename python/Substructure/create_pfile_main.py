import numpy as np
import snapshot
import h5py
import DistanceTool as distance
import matplotlib.pyplot as plt
import groupcat

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
#dm_ms =  header['MassTable'][1] # for Dandan
#dm_ms =  header['MassTable'][1]*1e10

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor
#boxsize = boxsize/1000./cosmopara.h # cMpc, glamer
boxsize = boxsize*a/1000./cosmopara.h # Mpc
#boxsize = boxsize# ckpc/h


## read in Galaxy catalog
catalog = '../../data/illustris_1/Galaxy_0'+str(snapNum)+'_test.dat'

GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
#CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])/1000./cosmopara.h # cMpc, glamer needs
#CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])/1000./cosmopara.h # cMpc
#CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])/1000./cosmopara.h # cMpc

CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])*a/1000./cosmopara.h # Mpc
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])*a/1000./cosmopara.h # Mpc
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])*a/1000./cosmopara.h # Mpc

#CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1]) # ckpc/h
#CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2]) # ckpc/h
#CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3]) # ckpc/h


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

subID = 218142

firstsubID=groupcat.loadHalos(basePath,snapNum, fields = ['GroupFirstSub'])
mainID=list(firstsubID).index(subID)

idx = list(GalaxyID).index(subID)
print idx

halo_dm = snapshot.loadHalo(basePath,snapNum,mainID,'dm')
#coord = halo_dm['Coordinates']/1000./cosmopara.h # cMpc, glamer
coord = halo_dm['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_dm['Coordinates']
dm_x,dm_y,dm_z = coord[:,0],coord[:,1],coord[:,2]

if (max(dm_x)-min(dm_x)> 0.5*boxsize): 
	dm_x = boundary(dm_x)

if (max(dm_y)-min(dm_y)> 0.5*boxsize): 
	dm_y = boundary(dm_y)

if (max(dm_z)-min(dm_z)> 0.5*boxsize): 
	dm_z = boundary(dm_z)

print 'halo dm coord loaded'

halo_gas = snapshot.loadHalo(basePath,snapNum,mainID,'gas')
#coord = halo_gas['Coordinates']/1000./cosmopara.h # cMpc,glamer
coord = halo_gas['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_gas['Coordinates']
gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
	gas_x = boundary(gas_x)
if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
	gas_y = boundary(gas_y)
if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
	gas_z = boundary(gas_z)

print 'halo gas coord loaded'

halo_st = snapshot.loadHalo(basePath,snapNum,mainID,'stars')
#coord = halo_st['Coordinates']/1000./cosmopara.h # cMpc,glamer
coord = halo_st['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_st['Coordinates']
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

if (max(st_x)-min(st_x)> 0.5*boxsize): 
	st_x = boundary(st_x)
if (max(st_y)-min(st_y)> 0.5*boxsize): 
	st_y = boundary(st_y)
if (max(st_z)-min(st_z)> 0.5*boxsize): 
	st_z = boundary(st_z)	

print 'halo stars coord loaded'


## ---- subhalo particles ------- ##

subhalo_dm = snapshot.loadSubhalo(basePath,snapNum,subID,'dm')
#coord = subhalo_dm['Coordinates']/1000./cosmopara.h # cMpc, glamer
coord = subhalo_dm['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_dm['Coordinates']
sub_dm_x,sub_dm_y,sub_dm_z = coord[:,0],coord[:,1],coord[:,2]

if (max(sub_dm_x)-min(sub_dm_x)> 0.5*boxsize): 
	sub_dm_x = boundary(sub_dm_x)

if (max(sub_dm_y)-min(sub_dm_y)> 0.5*boxsize): 
	sub_dm_y = boundary(sub_dm_y)

if (max(sub_dm_z)-min(sub_dm_z)> 0.5*boxsize): 
	sub_dm_z = boundary(sub_dm_z)

print 'subhalo dm coord loaded'

subhalo_gas = snapshot.loadSubhalo(basePath,snapNum,subID,'gas')
#coord = subhalo_gas['Coordinates']/1000./cosmopara.h # cMpc,glamer
coord = subhalo_gas['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_gas['Coordinates']
sub_gas_x,sub_gas_y,sub_gas_z = coord[:,0],coord[:,1],coord[:,2]

if (max(sub_gas_x)-min(sub_gas_x)> 0.5*boxsize): 
	sub_gas_x = boundary(sub_gas_x)
if (max(sub_gas_y)-min(sub_gas_y)> 0.5*boxsize): 
	sub_gas_y = boundary(sub_gas_y)
if (max(sub_gas_z)-min(sub_gas_z)> 0.5*boxsize): 
	sub_gas_z = boundary(sub_gas_z)

print 'subhalo gas coord loaded'

subhalo_st = snapshot.loadSubhalo(basePath,snapNum,subID,'stars')
#coord = subhalo_st['Coordinates']/1000./cosmopara.h # cMpc,glamer
coord = subhalo_st['Coordinates']*a/1000./cosmopara.h # Mpc
#coord = subhalo_st['Coordinates']
sub_st_x,sub_st_y,sub_st_z = coord[:,0],coord[:,1],coord[:,2]

if (max(sub_st_x)-min(sub_st_x)> 0.5*boxsize): 
	sub_st_x = boundary(sub_st_x)
if (max(sub_st_y)-min(sub_st_y)> 0.5*boxsize): 
	sub_st_y = boundary(sub_st_y)
if (max(sub_st_z)-min(sub_st_z)> 0.5*boxsize): 
	sub_st_z = boundary(sub_st_z)	

print 'subhalo stars coord loaded'


## change ref point to subhalo CM

dm_x,dm_y,dm_z = dm_x-CM_x[idx],dm_y-CM_y[idx],dm_z-CM_z[idx]
gas_x,gas_y,gas_z = gas_x-CM_x[idx],gas_y-CM_y[idx],gas_z-CM_z[idx]
st_x,st_y,st_z = st_x-CM_x[idx],st_y-CM_y[idx],st_z-CM_z[idx]

sub_dm_x,sub_dm_y,sub_dm_z = sub_dm_x-CM_x[idx],sub_dm_y-CM_y[idx],sub_dm_z-CM_z[idx]
sub_gas_x,sub_gas_y,sub_gas_z = sub_gas_x-CM_x[idx],sub_gas_y-CM_y[idx],sub_gas_z-CM_z[idx]
sub_st_x,sub_st_y,sub_st_z = sub_st_x-CM_x[idx],sub_st_y-CM_y[idx],sub_st_z-CM_z[idx]
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

sub_dm_x,sub_dm_y,sub_dm_z = sub_dm_x-ref_x,sub_dm_y-ref_y,sub_dm_z-ref_z
sub_gas_x,sub_gas_y,sub_gas_z = sub_gas_x-ref_x,sub_gas_y-ref_y,sub_gas_z-ref_z
sub_st_x,sub_st_y,sub_st_z = sub_st_x-ref_x,sub_st_y-ref_y,sub_st_z-ref_z


# mass
gas_ms = halo_gas['Masses']*1e10/cosmopara.h #Msun, glamer
st_ms = halo_st['Masses']*1e10/cosmopara.h
sub_gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h #Msun, glamer
sub_st_ms = subhalo_st['Masses']*1e10/cosmopara.h
#gas_ms = subhalo_gas['Masses'] #10^10 Msun/h
#st_ms = subhalo_st['Masses']

## ------ seperate substructure particles out ------ ##

mask_dm=np.in1d(dm_x,sub_dm_x,invert=True)
mask_st=np.in1d(st_x,sub_st_x,invert=True)
mask_gas=np.in1d(gas_x,sub_gas_x,invert=True)

osub_dm_x,osub_dm_y,osub_dm_z=dm_x[mask_dm],dm_y[mask_dm],dm_z[mask_dm]
osub_gas_x,osub_gas_y,osub_gas_z=gas_x[mask_gas],gas_y[mask_gas],gas_z[mask_gas]
osub_st_x,osub_st_y,osub_st_z=st_x[mask_st],st_y[mask_st],st_z[mask_st]

osub_gas_ms=gas_ms[mask_gas]
osub_st_ms=st_ms[mask_st]

### write in mass file 

#out_file = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'osub_dm.dat','w')
#out_file2 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'osub_gas.dat','w')
#out_file3 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'osub_st.dat','w')

out_file = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'wsub_dm_test.dat','w')
out_file2 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'wsub_gas_test.dat','w')
out_file3 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'wsub_st_test.dat','w')

out_file4 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'_dm_test.dat','w')
out_file5 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'_gas_test.dat','w')
out_file6 = open('../../data/illustris_1/snapshot_'+str(snapNum)+'_particle/particle_'+str(subID)+'_st_test.dat','w')

out_file.write('# nparticles '+str(dm_x.size)+'\n')
out_file2.write('# nparticles '+str(gas_x.size)+'\n')
out_file3.write('# nparticles '+str(st_x.size)+'\n')
out_file4.write('# nparticles '+str(sub_dm_x.size)+'\n')
out_file5.write('# nparticles '+str(sub_gas_x.size)+'\n')
out_file6.write('# nparticles '+str(sub_st_x.size)+'\n')
'''
for i in range(osub_dm_x.size):
	out_file.write(str(osub_dm_x[i])+'\t'+str(osub_dm_y[i])+'\t'+str(osub_dm_z[i])+'\t'+str(dm_ms)+'\n')

for i in range(osub_gas_x.size):
	out_file2.write(str(osub_gas_x[i])+'\t'+str(osub_gas_y[i])+'\t'+str(osub_gas_z[i])+'\t'+str(osub_gas_ms[i])+'\n')

for i in range(osub_st_x.size):
	out_file3.write(str(osub_st_x[i])+'\t'+str(osub_st_y[i])+'\t'+str(osub_st_z[i])+'\t'+str(osub_st_ms[i])+'\n')
'''

for i in range(dm_x.size):
	out_file.write(str(dm_x[i])+'\t'+str(dm_y[i])+'\t'+str(dm_z[i])+'\t'+str(dm_ms)+'\n')

for i in range(gas_x.size):
	out_file2.write(str(gas_x[i])+'\t'+str(gas_y[i])+'\t'+str(gas_z[i])+'\t'+str(gas_ms[i])+'\n')

for i in range(st_x.size):
	out_file3.write(str(st_x[i])+'\t'+str(st_y[i])+'\t'+str(st_z[i])+'\t'+str(st_ms[i])+'\n')

for i in range(sub_dm_x.size):
	out_file4.write(str(sub_dm_x[i])+'\t'+str(sub_dm_y[i])+'\t'+str(sub_dm_z[i])+'\t'+str(dm_ms)+'\n')

for i in range(sub_gas_x.size):
	out_file5.write(str(sub_gas_x[i])+'\t'+str(sub_gas_y[i])+'\t'+str(sub_gas_z[i])+'\t'+str(sub_gas_ms[i])+'\n')

for i in range(sub_st_x.size):
	out_file6.write(str(sub_st_x[i])+'\t'+str(sub_st_y[i])+'\t'+str(sub_st_z[i])+'\t'+str(sub_st_ms[i])+'\n')
