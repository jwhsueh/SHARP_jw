import numpy as np
import snapshot
import h5py
import DistanceTool as distance
import matplotlib.pyplot as plt

import pandas as pd

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1/'
snapNum = 99

class cosmopara:
	h = 0.704
	OM = 0.27

## read in Galaxy catalog

catalog_x = '../../data/illustris_1/full_'+str(snapNum)+'_x.dat'
catalog_y = '../../data/illustris_1/full_'+str(snapNum)+'_y.dat'
catalog_z = '../../data/illustris_1/full_'+str(snapNum)+'_z.dat'

table_x,table_y,table_z = pd.read_csv(catalog_x,sep = '\t'),pd.read_csv(catalog_y,sep = '\t'),pd.read_csv(catalog_z,sep = '\t')
# add projection flag
table_x['proj'],table_y['proj'],table_z['proj'] = 0,1,2

## to a long table
table = table_x.append(table_y,ignore_index = True)
table = table.append(table_z,ignore_index = True)

## velocity dispersion pick
table = table.loc[(table.velDisp > 130) & (table.velDisp < 370),:]

## ------ snapshot info ----- ##

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

dm_ms =  header['MassTable'][1]*1e10/cosmopara.h # mass of dm particle [M_sun]
boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

class lenspara:
	zl = redshift
	zs = 2.0

# critical density
sigma_c = distance.critical_density(cosmopara,lenspara)*cosmopara.h # M_sun/Mpc^2
sigma_c = sigma_c/(distance.mpc2arcs(cosmopara,1.,lenspara.zl))**2 # M_sun/arcsec^2

## ----- simulation box functions ----- ##

def boundary(ci):
	
	#ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	for i in ci:
		if i<0.5*boxsize:
			i = i+boxsize

	return ci

## ---- convergence plot ----- ## 

GalaxyID = 31645
projection = 0

galaxy_data = table.loc[(table.subfindID == GalaxyID) & (table.proj == projection),:]
CM_x,CM_y,CM_z = float(galaxy_data.pos_x),float(galaxy_data.pos_y),float(galaxy_data.pos_z)

## create particle catalog
particle_dm,particle_gas,particle_st = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

subhalo_dm = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'dm')
coord = subhalo_dm['Coordinates'] # ckpc/h
dm_x,dm_y,dm_z = coord[:,0],coord[:,1],coord[:,2]

if (max(dm_x)-min(dm_x)> 0.5*boxsize): 
	dm_x = boundary(dm_x)

if (max(dm_y)-min(dm_y)> 0.5*boxsize): 
	dm_y = boundary(dm_y)

if (max(dm_z)-min(dm_z)> 0.5*boxsize): 
	dm_z = boundary(dm_z)


dm_x,dm_y,dm_z = dm_x - CM_x, dm_y - CM_y, dm_z - CM_z # ckpc/h
dm_x,dm_y,dm_z = dm_x*a/cosmopara.h,dm_y*a/cosmopara.h,dm_z*a/cosmopara.h # kpc
dm_x = distance.mpc2arcs(cosmopara,dm_x/1000.,lenspara.zl)
dm_y = distance.mpc2arcs(cosmopara,dm_y/1000.,lenspara.zl)
dm_z = distance.mpc2arcs(cosmopara,dm_z/1000.,lenspara.zl)

particle_dm['pos_x'],particle_dm['pos_y'],particle_dm['pos_z'] = dm_x,dm_y,dm_z # arcsec
particle_dm['mass'] = dm_ms

print 'dm particle loaded'

subhalo_gas = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'gas')
coord = subhalo_gas['Coordinates']
gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h

if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
	gas_x = boundary(gas_x)
if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
	gas_y = boundary(gas_y)
if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
	gas_z = boundary(gas_z)

gas_x,gas_y,gas_z = gas_x - CM_x, gas_y - CM_y, gas_z - CM_z # ckpc/h
gas_x,gas_y,gas_z = gas_x*a/cosmopara.h, gas_y*a/cosmopara.h, gas_z*a/cosmopara.h # kpc
gas_x = distance.mpc2arcs(cosmopara,gas_x/1000.,lenspara.zl)
gas_y = distance.mpc2arcs(cosmopara,gas_y/1000.,lenspara.zl)
gas_z = distance.mpc2arcs(cosmopara,gas_z/1000.,lenspara.zl)

particle_gas['pos_x'],particle_gas['pos_y'],particle_gas['pos_z'] = gas_x,gas_y,gas_z # arcsec
particle_gas['mass'] = gas_ms

print 'gas particle loaded'

subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'stars')
coord = subhalo_st['Coordinates']
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

st_ms = subhalo_st['Masses']*1e10/cosmopara.h

if (max(st_x)-min(st_x)> 0.5*boxsize): 
	gas_x = boundary(st_x)
if (max(st_y)-min(st_y)> 0.5*boxsize): 
	st_y = boundary(st_y)
if (max(st_z)-min(st_z)> 0.5*boxsize): 
	st_z = boundary(st_z)	

st_x,st_y,st_z = st_x - CM_x, st_y - CM_y, st_z - CM_z # ckpc/h
st_x,st_y,st_z = st_x*a/cosmopara.h, st_y*a/cosmopara.h, st_z*a/cosmopara.h # kpc
st_x = distance.mpc2arcs(cosmopara,st_x/1000.,lenspara.zl)
st_y = distance.mpc2arcs(cosmopara,st_y/1000.,lenspara.zl)
st_z = distance.mpc2arcs(cosmopara,st_z/1000.,lenspara.zl)

particle_st['pos_x'],particle_st['pos_y'],particle_st['pos_z'] = st_x,st_y,st_z # arcsec
particle_st['mass'] = st_ms

print 'stars particle loaded'

particle = particle_dm.append(particle_gas)
particle = particle.append(particle_st)

## ---- particle loading done --- ##

## create annulus catalog

r_cut = 3. # arcsec
delta = 0.0001
r_edge = np.linspace(0,r_cut,r_cut/delta+1)
r_edge = r_edge[1:]

annulus = pd.DataFrame()
annulus['r'] = r_edge - delta/2.

# projection distance
particle['rx_proj'] = np.sqrt(particle.pos_y**2+particle.pos_z**2)
particle['ry_proj'] = np.sqrt(particle.pos_x**2+particle.pos_z**2)
particle['rz_proj'] = np.sqrt(particle.pos_x**2+particle.pos_y**2)

tot_ms = np.zeros(r_edge.size)
tot_area = np.zeros(r_edge.size)

for i in range(r_edge.size):
	print i
	r_in = r_edge[i]
	inside = particle.loc[(particle.rx_proj<r_in),:]
	tot_ms[i] = np.sum(inside.mass)
	tot_area[i] = np.pi*r_in**2 

as_ms = np.zeros(r_edge.size)
as_area = np.zeros(r_edge.size)

for i in range(r_edge.size-1):
	as_ms[i+1] = tot_ms[i+1]-tot_ms[i]
	as_area[i+1] = tot_area[i+1]-tot_area[i]

annulus['mass'],annulus['area'] = as_ms,as_area
annulus['kappa'] = annulus.mass/annulus.area/sigma_c


# interp (x)

Einstein_radius = np.interp(0.5,annulus.r,annulus.kappa)
'''
i = 1
while True:
	diff1 = annulus.kappa[i-1]-0.6
	diif2 = annulus.kappa[i] - 0.6

	if diff1*diif2<0:
		print annulus.r[i-1],annulus.r[i]
		print annulus.kappa[i-1],annulus.kappa[i]
		Einstein_radius = (annulus.r[i-1]+annulus.r[i])/2.
		break
	else:
		i = i+1
'''

#r_peak = annulus.r.loc[annulus.kappa == np.max(annulus.kappa)]
#Einstein_radius = Einstein_radius -r_peak
print Einstein_radius
print galaxy_data.R_E

#print annulus.kappa[220:240] 

#plt.plot(r_edge,annulus.kappa)
#plt.show()