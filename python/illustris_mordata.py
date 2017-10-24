import numpy as np
import snapshot
import h5py
import groupcat
import matplotlib.pyplot as plt

#table2 = np.loadtxt('../data/galaxy.txt',skiprows=1)

table = np.loadtxt('../data/Illustris_1/AllGalaxy_099_x.dat',skiprows=1)
mor = table[:,10]
mass = table[:,4]
ID = table[:,0]
bf = table[:,2]
df= table[:,1]
ser = table[:,8]

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'
ssNum = 99

color1 = np.zeros(ID.size)
color2 = np.zeros(ID.size)
g = np.zeros(ID.size)

mag_u = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloStellarPhotometrics'])[:,1]
mag_g = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloStellarPhotometrics'])[:,4]
mag_r = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloStellarPhotometrics'])[:,3]

print ID.size

for i in range(ID.size):

	#if mag_g[ID[i]]<(-22):

	color1[i]=mag_u[ID[i]]-mag_g[ID[i]]
	color2[i]=mag_g[ID[i]]-mag_r[ID[i]]
	g[i] = mag_r[ID[i]]

np.savetxt('galaxy_table.txt',np.c_[ser,g,color2,mor])

'''
#plt.scatter(np.log10(mass[np.logical_and(mor==1,ser>2)]),color2[np.logical_and(mor==1,ser>2)],marker='.',color='b')
#plt.scatter(np.log10(mass[np.logical_and(mor==0,ser<2)]),color2[np.logical_and(mor==0,ser<2)],marker='.',color='r')

#plt.scatter(np.log10(mass[np.logical_and(bf<0.6,df>0.4)]),color2[np.logical_and(bf<0.6,df>0.4)],marker='.',color='b')
#plt.scatter(np.log10(mass[np.logical_and(bf>0.6,df<0.4)]),color2[np.logical_and(bf>0.6,df<0.4)],marker='.',color='r')

#plt.scatter(ser[np.logical_and(mor==1,ser>2)],color2[np.logical_and(mor==1,ser>2)],marker='.',color='b')
#plt.scatter(ser[np.logical_and(mor==0,ser<2)],color2[np.logical_and(mor==0,ser<2)],marker='.',color='r')

plt.scatter(g[mor==1],color2[mor==1],marker='.',color='b',label='spiral')
plt.scatter(g[mor==0],color2[mor==0],marker='.',color='r',label='elliptical')
plt.xlabel('abs magnitude')
plt.ylabel('color')
plt.title('galaxy morphology bimodial')
plt.legend()
#plt.show()
plt.savefig('galaxy_bimodial2.png')
'''
