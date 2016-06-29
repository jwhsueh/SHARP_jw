import numpy as np
import snapshot
import h5py
import DistanceTool as distance
import matplotlib.pyplot as plt
import scipy.ndimage

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

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

# region
r_cut = 30. # 30 kpc/h

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
	else: # x-y plane
		p0,p1 = c0,c1

	return p0,p1


## -------start to plot particles

ID = 93295

idx = list(GalaxyID).index(ID)
print idx
GalaxyID = ID

# Use star particles only

subhalo_st = snapshot.loadSubhalo(basePath,snapNum,GalaxyID,'stars')
coord = subhalo_st['Coordinates']
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

if (max(st_x)-min(st_x)> 0.5*boxsize): 
	gas_x = boundary(st_x)
if (max(st_y)-min(st_y)> 0.5*boxsize): 
	st_y = boundary(st_y)
if (max(st_z)-min(st_z)> 0.5*boxsize): 
	st_z = boundary(st_z)	

# change ref point to subhalo CM
st_x,st_y,st_z = st_x-CM_x[idx],st_y-CM_y[idx],st_z-CM_z[idx]

# distance mask
st_r = np.sqrt(st_x**2+st_y**2+st_z**2)

mask = st_r < r_cut

st_x,st_y,st_z = st_x[mask],st_y[mask],st_z[mask]

#st_x,st_y,st_z = a*st_x,a*st_y,a*st_z

# star particle masses
st_ms = subhalo_st['Masses']*1e10 # Msun/h
st_ms = st_ms[mask]

## grid & grid size

proj = 2 # projected axis

# projected

ax0,ax1 = projection(st_x,st_y,st_z,proj)

# assign to each grid point (0.1 kpc/h)
ax0,ax1 = np.round(ax0*10),np.round(ax1*10)

# get grid min & max
ax0_range = [np.min(ax0),np.max(ax0)]
ax1_range = [np.min(ax1),np.max(ax1)]

# create grid
grid_0 = np.arange(ax0_range[0],ax0_range[1]+1)
grid_1 = np.arange(ax1_range[0],ax1_range[1]+1)

mass_grid = np.zeros((grid_0.size,grid_1.size))

print grid_0.size,grid_1.size,mass_grid.shape
## decide which grid each particle falls into
## this take too much time
for i in range(ax0.size):
	print i
	idx_0 = list(grid_0).index(ax0[i])
	idx_1 = list(grid_1).index(ax1[i])

	mass_grid[idx_0,idx_1] = mass_grid[idx_0,idx_1]+st_ms[i]

# smoothing

mass_smooth = scipy.ndimage.filters.gaussian_filter(mass_grid,sigma = 3)
#mass_smooth = mass_smooth/np.max(mass_smooth)
#mass_smooth = np.log10(mass_smooth)

# plot 
#plt.scatter(st_x,st_z,marker='.')
CS = plt.contourf(grid_1/10,grid_0/10,mass_smooth,[2e6,4e6,6e6,1e7,1.4e7,1.8e7,2.2e7,2.6e7,3e7,5e7,7e7])
plt.colorbar()

plt.title('SubfindID 93295 , IncAngle = 28 deg')
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.xlabel('kpc/h')
plt.ylabel('kpc/h')
plt.show()