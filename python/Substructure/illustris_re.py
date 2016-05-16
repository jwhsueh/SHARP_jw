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
print dm_ms

## see if subhalo is on boundary

boxsize = header['BoxSize']  # ckpc/h
print boxsize

redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor
c = 3e5 # km/s

class lenspara:
	zl = redshift
	zs = 2.0

# lens critical density
sig_c = distance.critical_density(cosmopara,lenspara)*cosmopara.h # M_sun/Mpc^2
sig_c = sig_c/1e6 # M_sun/kpc^2
#print sig_c/cosmopara.h
#print distance.mpc2arcs(cosmopara,1.,redshift)

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

def root_find(f,arg,root,endpt,delta):
	s,e = endpt[0],endpt[1]
	m = (s+e)/2.

	while True:
		value1 = f(arg,s)
		value2 = f(arg,m)
		value3 = f(arg,e)

		print value1,value2,value3

		subs = np.array([value1-root,value2-root,value3-root])

		if np.abs(value2-root)<delta:
			bisec_root = m
			break

		elif subs[0]*subs[1]<0:
			e = m
			m = (s+e)/2.
		elif subs[1]*subs[2]<0:
			s = m
			m = (s+e)/2.
		else:
			print "root is outside searching range"
			bisec_root = np.NaN
			break


	return bisec_root

## surface mass density

def anulus_kappa(subhalo_class, r):
	step = 2 # kpc
	dm_dist,gas_dist,st_dist = subhalo_class.dm_d,subhalo_class.gas_d,subhalo_class.st_d

	dm_ms,gas_ms,st_ms = subhalo_class.dm_m,subhalo_class.gas_m,subhalo_class.st_m

	part_in,part_out = np.array([dm_dist<(r-step)]),np.array([dm_dist>(r+step)])
	mask = -(part_in+part_out)
	#mask = -(part_out)
	#print mask.size, dm_dist.size

	dm_part = dm_dist[list(mask)].size
	dm_tot = dm_part*dm_ms

	part_in,part_out = np.array([gas_dist<(r-step)]),np.array([gas_dist>(r+step)])
	mask = -(part_in+part_out)
	#mask = -(part_out)

	gas_part = gas_ms[list(mask)]
	gas_tot = np.sum(gas_part)

	part_in,part_out = np.array([st_dist<(r-step)]),np.array([st_dist>(r+step)])
	mask = -(part_in+part_out)
	#mask = -(part_out)

	st_part = st_ms[list(mask)]
	st_tot = np.sum(gas_part)

	tot_ms = dm_tot+gas_tot+st_tot
	#print tot_ms
	kappa = tot_ms/np.pi/((r+step)**2-(r-step)**2)/sig_c # M_sun/arcsec^2
	#kappa = tot_ms/np.pi/r**2/sig_c

	return kappa


############################

GalaxyID = np.array([105690,114445])
CM_x = np.array([14435.1,37873.5])
CM_y = np.array([60465.2,1752.11])
CM_z = np.array([71860.8,36323.1])

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

	#print min(np.abs(dm_x))

	## projection & R_e
	proj_axis = 0

	dm_p0,dm_p1 = projection(dm_x,dm_y,dm_z,proj_axis)
	gas_p0,gas_p1 = projection(gas_x,gas_y,gas_z,proj_axis)
	st_p0,st_p1 = projection(st_x,st_y,st_z,proj_axis)

	print min(np.abs(dm_p0))

	# distance to CM
	dm_dist = np.sqrt(dm_p0**2+dm_p1**2) # ckpc/h [comoving]
	gas_dist = np.sqrt(gas_p0**2+gas_p1**2)
	st_dist = np.sqrt(st_p0**2+st_p1**2)

	dm_dist = dm_dist/cosmopara.h*a # kpc
	gas_dist = gas_dist/cosmopara.h*a # kpc
	st_dist = st_dist/cosmopara.h*a # kpc

	print min(dm_dist)

	#dm_dist = distance.mpc2arcs(cosmopara,dm_dist,redshift) # arcsec
	#gas_dist = distance.mpc2arcs(cosmopara,gas_dist,redshift)
	#st_dist = distance.mpc2arcs(cosmopara,st_dist,redshift)

	# mass
	gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h
	st_ms = subhalo_st['Masses']*1e10/cosmopara.h

	class subhalo_info:
		dm_d,gas_d,st_d = dm_dist,gas_dist,st_dist
		dm_m,gas_m,st_m = dm_ms,gas_ms,st_ms

	# bisection to find R_e
	endpt = np.array([1,20])
	delta = 0.01

	r = np.linspace(0,20,201)
	kappa= []
	for ri in r:
		kappa.append(anulus_kappa(subhalo_info,ri))

	#r,step = 7,1
	#part_in,part_out = np.array([dm_dist<(r-step)]),np.array([dm_dist>(r+step)])
	#mask = -(part_in+part_out)

	#x,y = dm_p0[list(mask)],dm_p1[list(mask)]
	#plt.scatter(x,y)


	plt.plot(r,kappa)
	plt.xlabel('kpc')
	plt.ylabel('surface mass density')
	plt.show()

	#R_e = root_find(anulus_kappa,subhalo_info,1.0,endpt,delta) #Mpc
	#R_e = distance.mpc2arcs(cosmopara,R_e,redshift)
	#print R_e
			

	#R_e = np.interp(sig_c,sigma,)

	## no we shouldn't do interp
	## this is a root finding prob



