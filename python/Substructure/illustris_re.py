import numpy as np
import snapshot
import h5py
import DistanceTool as distance

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = 99

class cosmopara:
	h = 0.704
	OM = 0.27

c = 3e8 # m/s

## read in Galaxy catalog

catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])

p_type = ['gas','dm','stars']

field = ['Masses']

## read in dm particle mass

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

dm_ms =  header['MassTable'][1]*1e10/cosmopara.h # mass of dm particle [M_sun]

## see if subhalo is on boundary

boxsize = header['BoxSize']  # ckpc/h
print boxsize

redshift = header['Redshift']

class lenspara:
	zl = redshift
	zs = 2.0

# lens critical density
sig_c = distance.critical_density(cosmopara,lenspara) # M_sun/Mpc^2
sig_c = sig_c/(distance.mpc2arcs(cosmopara,1.0,redshift))**2.  # M_sun/arcsec^2

## this function deal w/ galaxy on the boundary
def boundary(ci):
	
	ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	return ci

## project along one axis
def projection(c0,c1,c2,axis):

	if axis == 0: # y-z plane
		p0,p1 = c1,c2
	elif axis == 1: # x-z plane
		p0,p1 = c0,c2
	else: # x-y plane
		p0,p1 = c0,c1

	return p0,p1

## this is for bisection method

def bisection(f,arg,root,endpt,delta):
	s,e = endpt[0],endpt[1]

	while True:
		value1 = f(arg,s)
		value2 = f(arg,e)

		if np.abs(value1-root)<delta:
			bisec_root = value1
			break

		elif np.abs(value2-root)<delta:
			bisec_root = value2
			break

		elif np.abs(value1-root)<np.abs(value2-root):
			e = (s+e)/2.0
		else:
			s = (s+e)/2.0




GalaxyID = GalaxyID[:2]

for i in range(GalaxyID.size):

	print GalaxyID[i]

	subhalo_dm = snapshot.loadSubhalo(basePath,snapNum,GalaxyID[i],'dm')
	coord = subhalo_dm['Coordinates'] # ckpc/h
	dm_x,dm_y,dm_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(dm_x)-min(dm_x)> 0.5*boxsize): 
		dm_x = boundary(dm_x)
	if (max(dm_y)-min(dm_y)> 0.5*boxsize): 
		dm_y = boundary(dm_y)
	if (max(dm_z)-min(dm_z)> 0.5*boxsize): 
		dm_z = boundary(dm_z)

	print 'dm coord loaded'

	subhalo_gas = snapshot.loadSubhalo(basePath,snapNum,GalaxyID[i],'gas')
	coord = subhalo_gas['Coordinates']
	gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
		gas_x = boundary(gas_x)
	if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
		gas_y = boundary(gas_y)
	if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
		gas_z = boundary(gas_z)

	print 'gas coord loaded'

	subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID[i],'stars')
	coord = subhalo_st['Coordinates']
	st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(st_x)-min(st_x)> 0.5*boxsize): 
		gas_x = boundary(st_x)
	if (max(st_y)-min(st_y)> 0.5*boxsize): 
		st_y = boundary(st_y)
	if (max(st_z)-min(st_z)> 0.5*boxsize): 
		st_z = boundary(st_z)	

	print 'stars coord loaded'

	## change ref point to subhalo CM

	dm_x,dm_y,dm_z = dm_x-CM_x[i],dm_y-CM_y[i],dm_z-CM_z[i]
	gas_x,gas_y,gas_z = gas_x-CM_x[i],gas_y-CM_y[i],gas_z-CM_z[i]
	st_x,st_y,st_z = st_x-CM_x[i],st_y-CM_y[i],st_z-CM_z[i]

	## projection & R_e
	proj_axis = 2

	dm_p0,dm_p1 = projection(dm_x,dm_y,dm_z,proj_axis)
	gas_p0,gas_p1 = projection(gas_x,gas_y,gas_z,proj_axis)
	st_p0,st_p1 = projection(st_x,st_y,st_z,proj_axis)

	# distance to CM
	dm_dist = np.sqrt(dm_p0**2+dm_p1**2) # ckpc/h
	gas_dist = np.sqrt(gas_p0**2+gas_p1**2)
	st_dist = np.sqrt(st_p0**2+st_p1**2)

	dm_dist = dm_dist*c/cosmopara.h/1000. # Mpc
	gas_dist = gas_dist*c/cosmopara.h/1000. # Mpc
	st_dist = st_dist*c/cosmopara.h/1000. # Mpc

	dm_dist = distance.mpc2arcs(cosmopara,dm_dist,redshift) # arcsec
	gas_dist = distance.mpc2arcs(cosmopara,gas_dist,redshift)
	st_dist = distance.mpc2arcs(cosmopara,st_dist,redshift)

	# mass
	gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h
	st_ms = subhalo_st['Masses']*1e10/cosmopara.h

	step = max(st_dist)*100. # resolution 0.01"
	dist_range = np.linspace(0,max(st_dist),step)

	sigma = np.zeros(dist_range.size)
	inter = 0

	for j in range(dist_range.size):
		# dm
		dm_part = dm_dist[dm_dist<dist_range[j]].size
		dm_tot = dm_part*dm_ms

		# gas
		gas_part = gas_ms[gas_dist<dist_range[j]]
		gas_tot = np.sum(gas_part)

		# stars
		st_part = st_ms[st_dist<dist_range[j]]
		st_tot = np.sum(gas_part)

		tot_ms = dm_tot+gas_tot+st_tot
		sigma[j] = tot_ms/np.pi/dist_range[j]**2. # M_sun/arcsec^2

		if sigma[j]>sig_c:
			inter = inter+1

		if inter == 10:
			break

	#R_e = np.interp(sig_c,sigma,)

	## no we shouldn't do interp
	## this is a root finding prob



