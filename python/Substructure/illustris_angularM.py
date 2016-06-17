import numpy as np
import snapshot
import h5py
import groupcat
import DistanceTool as distance

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = 99

class cosmopara:
	h = 0.704
	OM = 0.27

## read in Galaxy catalog

catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'_str.dat'
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

# get peculiar velocity of group
SubVel = groupcat.loadSubhalos(basePath,snapNum, fields = ['SubhaloVel'])

## this function deal w/ galaxy on the boundary
def boundary(ci):
	
	#ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	for i in ci:
		if i<0.5*boxsize:
			i = i+boxsize

	return ci

## calculate angular momentum on each axis

Axis = np.zeros((GalaxyID.size,3))*np.nan  # principle ratational axis of each galaxy
theta_i = np.zeros((GalaxyID.size,3))*np.nan # inclination angle of three main axis

in_cat = open('../../data/illustris_1/Inclination_'+str(snapNum)+'_str.dat','w')
in_cat.write('# Galaxy ID 	theta_x		theta_y		theta_z \n')

for i in range(GalaxyID.size):
	print GalaxyID[i]

	subID = GalaxyID[i] 
	GalaxyVel = SubVel[subID]

	# only use star particles

	subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID[i],'stars')
	coord = subhalo_st['Coordinates']
	st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(st_x)-min(st_x)> 0.5*boxsize): 
		gas_x = boundary(st_x)
	if (max(st_y)-min(st_y)> 0.5*boxsize): 
		st_y = boundary(st_y)
	if (max(st_z)-min(st_z)> 0.5*boxsize): 
		st_z = boundary(st_z)	

	vel = subhalo_st['Velocities']
	# relative velocity to center
	st_vx,st_vy,st_vz = vel[:,0]*np.sqrt(a)-GalaxyVel[0],vel[:,1]*np.sqrt(a)--GalaxyVel[1],vel[:,2]*np.sqrt(a)-GalaxyVel[2]

	st_ms = subhalo_st['Masses']*1e10/cosmopara.h

	st_Lx = np.sum((st_y*st_vz-st_z*st_vy)*st_ms)
	st_Ly = np.sum((st_z*st_vx-st_x*st_vz)*st_ms)
	st_Lz = np.sum((st_x*st_vy-st_y*st_vx)*st_ms)

	L_len = np.sqrt(st_Lx**2+st_Ly**2+st_Lz**2)

	L_xu,L_yu,L_zu = st_Lx/L_len,st_Ly/L_len,st_Lz/L_len  # unit vector component

	Axis[i,:] = np.array([L_xu,L_yu,L_zu])
	print
	theta_x,theta_y,theta_z = np.arccos(L_xu),np.arccos(L_yu),np.arccos(L_zu)
	theta_i[i,:] = np.array([np.degrees(theta_x),np.degrees(theta_y),np.degrees(theta_z)])
	print theta_i[i,:]

	in_cat.write(str(GalaxyID[i])+'	'+str(theta_i[i,:])[1:-1]+'\n')






