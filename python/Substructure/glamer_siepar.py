import numpy as np

filename = '/Volumes/sting_1/snap99_111/particle_227340_sie.dat'
table = np.loadtxt(filename)
#print len(table[:,0])
mass = np.empty(len(table[:,0]))
#print mass
mass.fill(6.26273e6)
#print mass

arr = np.empty((len(table[:,0]),4))
arr[:,:3] = table
arr[:,3] =mass

np.savetxt(filename,arr)