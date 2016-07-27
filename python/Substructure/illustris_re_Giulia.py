import numpy as np
import snapshot
import h5py
import DistanceTool as distance
import matplotlib.pyplot as plt

import pandas as pd

####### change path here #######
### basePath = snapshot dir
### catalog = Galaxy catalog

snapNum = 99
basePath = '/Volumes/narsil_1/jwhsueh/illustris_1/'
catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'_sig.dat'

output_name = basePath+'snapshot'+str(snapNum)+'_Re.dat'

###### ------don't need to change anything below----- ########

class cosmopara:
	h = 0.704
	OM = 0.27

## read in Galaxy catalog

#cata_field = ['subfindID','pos_x','pos_y','pos_z','mass','stellar_mass','velDisp']
cata_field = ['subfindID','pos_x','pos_y','pos_z','mass','stellar_mass','velDisp','7','8','9']
table = pd.read_csv(catalog,sep = '\s+',comment = '#',names = cata_field)


## velocity dispersion pick
#table = table.loc[(table.velDisp > 130) & (table.velDisp < 370),:]

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

## ---- calculating Einstein radius ----- ## 

subfind_ID = np.array(table.subfindID)
output_file = open(output_name,'w')
output_file.write('subfindID\tR_Ex\tR_Ey\tR_Ez\n')

Re_proj = np.zeros((len(subfind_ID),3))

k = 0
for GalaxyID in subfind_ID:
	print 'subfindID = '+str(GalaxyID)

	galaxy_data = table.loc[(table.subfindID == GalaxyID),:]
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

	print 'star particle loaded'

	particle = particle_dm.append(particle_gas)
	particle = particle.append(particle_st)

	## ---- particle loading done --- ##

	

	# projection distance
	particle['rx_proj'] = np.sqrt(particle.pos_y**2+particle.pos_z**2)
	particle['ry_proj'] = np.sqrt(particle.pos_x**2+particle.pos_z**2)
	particle['rz_proj'] = np.sqrt(particle.pos_x**2+particle.pos_y**2)

	## ---- calculate R_E for each projection ---- ##



	for projection in range(3):

		print 'projection = '+str(projection)


		# choose projection axis distance
		if projection == 0:
			proj_r = particle['rx_proj']
		elif projection == 1:
			proj_r = particle['ry_proj']
		else:
			proj_r = particle['rz_proj']

		print 'calculate Einstein radius...'

		rerange = 0.5

		## create annulus mock

		r_cut = 3. # arcsec
		delta = 0.01
		r_edge = np.linspace(0,r_cut,r_cut/delta+1)
		r_edge = r_edge[1:]

		#annulus = pd.DataFrame()
		as_r = r_edge - delta/2.

		## ---- iteration loop ---- ##
		for j in range(2):
			tot_ms = np.zeros(r_edge.size)
			tot_area = np.zeros(r_edge.size)

			for i in range(r_edge.size):
				r_in = r_edge[i]
				inside = particle.loc[(proj_r<r_in),:]
				tot_ms[i] = np.sum(inside.mass)
				tot_area[i] = np.pi*r_in**2 

			as_ms = np.zeros(r_edge.size)
			as_area = np.zeros(r_edge.size)

			for i in range(r_edge.size-1):
				as_ms[i+1] = tot_ms[i+1]-tot_ms[i]
				as_area[i+1] = tot_area[i+1]-tot_area[i]

			as_kappa = as_ms/as_area/sigma_c

			if j>0:
				break
		## ---------done iteration ----- ##

		Re = np.interp(0.5,as_r,as_kappa)
		## re-search R_E
		delta = delta/5.
		r_start,r_end = Re-rerange,Re+rerange

		if r_start<0:
			r_start = 0.

		if r_end >r_cut:
			r_end = r_cut

		r_edge = np.linspace(r_start,r_end,2.*rerange/delta+1)
		r_edge = r_edge[1:]
		as_r = r_edge - delta/2.

		for i in range(len(as_r)-1):
			if (as_kappa[i] - 0.5) <0:
				Re = as_r[i]
				break

		print 'R_Ein = '+str(Re)

		####------- end of calculation

		Re_proj[k,projection] = Re
		#print Re_proj[k,:]

		## write in the middle

	## Done w/ projections

	output_file.write(str(GalaxyID)+'\t'+str(Re_proj[k,:])[1:-1]+'\n')
	k = k+1

		
output_file.close()


