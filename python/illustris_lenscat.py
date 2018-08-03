import numpy as np
import h5py

basePath = '/Volumes/narsil_1/jwhsueh/illustris_1/'

catfile = h5py.File(basePath+'pmsd.hdf5','r')

s89 = catfile['Snapshot_89']
subid =  s89['SubfindID']
Re_x,Re_y,Re_z = s89['Rein_x'],s89['Rein_y'],s89['Rein_z']
#print Re_x[:10]

np.savetxt(basePath+'snapshot89_re.dat',np.c_[subid,Re_x,Re_y,Re_z],fmt='%d %f %f %f')