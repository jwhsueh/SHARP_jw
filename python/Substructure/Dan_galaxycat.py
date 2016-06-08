import numpy as np

# galaxy cat2

basePath = '../../data/illustris_1/'
ssNumber = 99

catalog = basePath+'Galaxy_'+str(ssNumber)+'.dat'
Galaxy_ID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

catalog2 = basePath+'Dandan_Lens'+str(ssNumber)+'.dat'
LenID = np.loadtxt(catalog2,dtype = 'int',unpack=True, usecols=[0])
Len_table = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[1:])


catalog3 = basePath+'Dandan_Lens'+str(ssNumber)+'.dat'
PhotID = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Phot_table = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[1:])

# create table
table = np.zeros((Len_table.shape[0]+Phot_table.shape[0]-2,len(Galaxy_ID)))

