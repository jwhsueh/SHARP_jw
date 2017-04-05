import numpy as np
import matplotlib.pyplot as plt

listname ='/Volumes/sting_1/data/snap99_tri2_Rcusp.txt'
listname2 = '/Volumes/sting_1/data/snap99_elp2_Rcusp.txt'
relist = '/Volumes/sting_1/data/snap99_tri2_size.txt'
relist2 = '/Volumes/sting_1/data/snap99_elp2_size.txt'
galaxy_file = '../../data/illustris_1/AllGalaxy_099_x.dat'
tab = np.loadtxt(galaxy_file,skiprows=1)
galaxy_id = tab[:,0]
mstar=tab[:,6]

tab = np.loadtxt(relist)
re = tab[:,1]

id_list = np.genfromtxt(listname,dtype='str')

ms_list = []
for lensname in id_list:
	subid = int(lensname.split('_')[0])
	print subid
	idx = list(galaxy_id).index(subid)
	ms_list.append(mstar[idx])

ms_list = np.array(ms_list)

mask = re>0.3
ms_list,re = ms_list[mask],re[mask]

tab = np.loadtxt(relist2)
re2 = tab[:,1]

id_list2 = np.genfromtxt(listname2,dtype='str')

ms_list2 = []
for lensname in id_list2:
	subid = int(lensname.split('_')[0])
	print subid
	idx = list(galaxy_id).index(subid)
	ms_list2.append(mstar[idx])

ms_list2 = np.array(ms_list2)

plt.scatter(np.log10(ms_list2),re2,color = 'b',label='elliptical')
plt.scatter(np.log10(ms_list),re,color = 'r',label='disc')
plt.legend(loc=2,scatterpoints=1)
plt.xlabel('log(M*)')
plt.ylabel('Einstein radius (")')
#plt.gca().set_aspect('equal')
plt.savefig('../../data/glamer/glamer_sm_scatter.png')
#plt.show()