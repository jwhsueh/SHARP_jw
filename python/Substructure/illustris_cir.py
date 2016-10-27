import numpy as np
import snapshot
import h5py
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt

'''
This code calculate circularity of each star particle

J(E) curve is determined by gasous particles

E: binding energy, E=U+K
'''

#### --------------change path setting here -------------#####
## snapshotPath: where the snapshot directory is
## catalog: Galaxy catalog (read in)
## in_cata: Inclination catalog (write out)

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'
ssNum = 99

catalog = '../../data/illustris_1/Galaxy_'+str(snapNum)+'_sig.dat'
catalog_th = '../../data/illustris_1/Inclination_'+str(snapNum)+'_HD.dat'

#### -------------------------------- #####

class cosmopara:
	h = 0.704
	OM = 0.27

## this function deal w/ galaxy on the boundary
def boundary(ci):
	
	#ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	for i in ci:
		if i<0.5*boxsize:
			i = i+boxsize

	return ci

##### ------------------------------- #####

## read in Galaxy catalog

GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
CM_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
CM_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])
R_half = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[7]) # kpc/h

p_type = ['gas','dm','stars']

field = ['Masses','Velocities']

## read in snapshot header

f = h5py.File(snapshot.snapPath(snapshotPath,ssNum),'r')
header = dict( f['Header'].attrs.items() )

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

# get peculiar velocity of group
SubVel = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloVel']) # km/s

## ----------- gas particles

subfindID = 72224

idx = list(GalaxyID).index(subfindID) 

sub_gas = snapshot.loadSubhalo(snapshotPath,snapNum,subfindID,'gas')

coord = sub_gas['Coordinates']
gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

if (max(gas_x)-min(gas_x)> 0.5*boxsize): 
	gas_x = boundary(gas_x)
if (max(gas_y)-min(gas_y)> 0.5*boxsize): 
	gas_y = boundary(gas_y)
if (max(gas_z)-min(gas_z)> 0.5*boxsize): 
	gas_z = boundary(gas_z)

## change ref point to subhalo CM
gas_x,gas_y,gas_z = gas_x-CM_x[idx],gas_y-CM_y[idx],gas_z-CM_z[idx]	# ckpc/h
gas_x,gas_y,gas_z = gas_x*a,gas_y*a,gas_z*a # kpc/h

## Hubble drag
hd_x,hd_y,hd_z = gas_x*100.*distance.Ez(cosmopara,redshift),gas_y*100.*distance.Ez(cosmopara,redshift),gas_z*100.*distance.Ez(cosmopara,redshift)

vel = sub_gas['Velocities']

gas_vx,gas_vy,gas_vz = vel[:,0]*np.sqrt(a),vel[:,1]*np.sqrt(a),vel[:,2]*np.sqrt(a) # km/s
gas_vx,gas_vy,gas_vz = gas_vx+hd_x,gas_vy+hd_y,gas_vz+hd_z # apply Hubble drag

## ----- velocity is done



#gas_ms = sub_gas['Masses']*1e10/cosmopara.h
#gas_ms = sub_gas['Masses']/cosmopara.h

## ----- kinematics energy
KE = 0.5*(gas_vz**2+gas_vy**2+gas_vz**2)

## ----- potential energy
U = sub_gas['Potential']/a # (km/s)^2
#U = sub_gas['Potential']/a*gas_ms # (km/s)^2

## ---- mask out particles in outer region

## only use particles w/i 2x half mass radius
gas_r = np.sqrt(gas_x**2+gas_y**2+gas_z**2)
r_cut = R_half[idx]/cosmopara.h
	
if r_cut == 0.0:
	r_cut = 10.

mask = gas_r < 2.0*r_cut

gas_x,gas_y,gas_z = gas_x[mask],gas_y[mask],gas_z[mask]
gas_vx,gas_vy,gas_vz = gas_vx[mask],gas_vy[mask],gas_vz[mask]
U,KE=U[mask],KE[mask]

## ----- binding energy
E = U+KE
#E=U

## ---------- principal axis anugular momentum ------------- ##

th_x = np.loadtxt(catalog_th,dtype = 'float',unpack=True, usecols=[1])
th_y = np.loadtxt(catalog_th,dtype = 'float',unpack=True, usecols=[1])
th_z = np.loadtxt(catalog_th,dtype = 'float',unpack=True, usecols=[1])

th_x,th_y,th_z = np.radians(th_x[idx]),np.radians(th_y[idx]),np.radians(th_z[idx])

## angular momtentum for each gas particle
#gas_Lx = (gas_y*gas_vz-gas_z*gas_vy)*gas_ms
#gas_Ly = (gas_z*gas_vx-gas_x*gas_vz)*gas_ms
#gas_Lz = (gas_x*gas_vy-gas_y*gas_vx)*gas_ms

gas_Lx = (gas_y*gas_vz-gas_z*gas_vy)
gas_Ly = (gas_z*gas_vx-gas_x*gas_vz)
gas_Lz = (gas_x*gas_vy-gas_y*gas_vx)

## project to pricipal axis
gas_Jz = gas_Lx*np.cos(th_x)+gas_Ly*np.cos(th_y)+gas_Lz*np.cos(th_z)

#plt.scatter(KE/1e5,gas_Jz/1e3)
plt.scatter(U/2/1e5,gas_Jz/1e3,marker='.',s=1)
plt.xlabel("Binding Energy E/(10^5 km^2/s^2)")
plt.ylabel("Jz/(10^3 kpc km/s)")
plt.title('Gas Particles')
plt.show()


