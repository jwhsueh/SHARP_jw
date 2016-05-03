import numpy as np
import snapshot
import h5py

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = 99

h = 0.704


## read in Galaxy catalog

catalog = '../../data/illustris_1/Galaxy_0'+str(snapNum)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
CM_x = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[1])
CM_y = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[2])
CM_z = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[3])

p_type = ['gas','dm','stars']

field = ['Masses']

## read in dm particle mass

f = h5py.File(snapshot.snapPath(basePath,snapNum),'r')
header = dict( f['Header'].attrs.items() )

dm_ms =  header['MassTable'][1]*1e10/h # mass of dm particle [M_sun]



