import numpy as np
import matplotlib.pyplot as plt
'''
all_list1 = np.loadtxt('../../data/illustris_1/snap99_p1_all.txt')
all_list2 = np.loadtxt('../../data/illustris_1/snap99_p2_all.txt')
all_list3 = np.loadtxt('../../data/illustris_1/snap99_p3_all.txt')


plt.scatter(np.log10(all_list1[:,0]),all_list1[:,1],marker='x',color='c',alpha=0.3,label='all lens')
plt.scatter(np.log10(all_list2[:,0]),all_list2[:,1],marker='x',color='c',alpha=0.3)
plt.scatter(np.log10(all_list3[:,0]),all_list3[:,1],marker='x',color='c',alpha=0.3)

# ellptical
e_list1 = np.loadtxt('../../data/illustris_1/snap99_p1_e.txt')
e_list2 = np.loadtxt('../../data/illustris_1/snap99_p2_e.txt')
e_list3 = np.loadtxt('../../data/illustris_1/snap99_p3_e.txt')
plt.scatter(np.log10(e_list1[:,0]),e_list1[:,1],marker='o',color='b',label='elliptical')
plt.scatter(np.log10(e_list2[:,0]),e_list2[:,1],marker='o',color='b')
plt.scatter(np.log10(e_list3[:,0]),e_list3[:,1],marker='o',color='b')
'''
# disc
d_list1 = np.loadtxt('../../data/illustris_1/snap99_p1_tri.txt')
d_list2 = np.loadtxt('../../data/illustris_1/snap99_p2_tri.txt')
d_list3 = np.loadtxt('../../data/illustris_1/snap99_p3_tri.txt')
#plt.scatter(np.log10(d_list1[:,3]),d_list1[:,2],marker='*',color='r',label='disc')
#plt.scatter(np.log10(d_list2[:,3]),d_list2[:,2],marker='*',color='r')
#plt.scatter(np.log10(d_list3[:,3]),d_list3[:,2],marker='*',color='r')
d_id1,d_id2,d_id3 = d_list1[:,0].astype(int),d_list2[:,0].astype(int),d_list3[:,0].astype(int)
d_re1,d_re2,d_re3 = d_list1[:,2],d_list2[:,2],d_list3[:,2]
d_dim1,d_dim2,d_dim3 = d_list1[:,3],d_list2[:,3],d_list3[:,3]

# get inclination
tab = np.loadtxt('../../data/illustris_1/inclination_099_sig.dat')
in_ID = tab[:,0].astype(int)
in1,in2,in3 = tab[:,1],tab[:,2],tab[:,3]

d_in1,d_in2,d_in3 = np.zeros(len(d_id1)),np.zeros(len(d_id2)),np.zeros(len(d_id3))
for i in range(len(d_id1)):
	if np.in1d(d_id1[i],in_ID):
		idx = list(in_ID).index(d_id1[i])
		d_in1[i] = in1[idx]

for i in range(len(d_id2)):
	if np.in1d(d_id2[i],in_ID):
		idx = list(in_ID).index(d_id2[i])
		d_in2[i] = in2[idx]

for i in range(len(d_id3)):
	if np.in1d(d_id3[i],in_ID):
		idx = list(in_ID).index(d_id3[i])
		d_in3[i] = in3[idx]

plt.scatter(np.log10(d_dim1),d_re1,marker='x',color='k',alpha=0.3,label='disc galaxy')
plt.scatter(np.log10(d_dim2),d_re2,marker='x',color='k',alpha=0.3)
plt.scatter(np.log10(d_dim2),d_re2,marker='x',color='k',alpha=0.3)

## proj1
f_mask=(d_in1<=50.)|(d_in1>=130.)
e_mask=(d_in1>=80.)&(d_in1<=100.)

plt.scatter(np.log10(d_dim1[f_mask]),d_re1[f_mask],marker='*',color='g',label='face-on')
plt.scatter(np.log10(d_dim1[e_mask]),d_re1[e_mask],marker='o',color='r',label='edge-on')

## proj2
f_mask=(d_in2<=50.)|(d_in2>=130.)
e_mask=(d_in2>=80.)&(d_in2<=100.)

plt.scatter(np.log10(d_dim2[f_mask]),d_re2[f_mask],marker='*',color='g')
plt.scatter(np.log10(d_dim2[e_mask]),d_re2[e_mask],marker='o',color='r')

## proj3
f_mask=(d_in3<=50.)|(d_in3>=130.)
e_mask=(d_in3>=80.)&(d_in3<=100.)

plt.scatter(np.log10(d_dim3[f_mask]),d_re3[f_mask],marker='*',color='g')
plt.scatter(np.log10(d_dim3[e_mask]),d_re3[e_mask],marker='o',color='r')



#plt.scatter(np.log10(ms_list2),re2,color = 'b',marker='^',label='elliptical')
#plt.scatter(np.log10(ms_list),re,color = 'r',marker='v',label='disc')
plt.legend(loc=2,scatterpoints=1)
plt.xlabel('log(M_disc)')
plt.ylabel('Einstein radius (")')
plt.ylim(0.1,1.5)
plt.xlim(9.6,11.1)
plt.plot([9,11.5],[0.3,0.3],linestyle='--',color='k',lw=2,alpha=0.7)
#plt.gca().set_aspect('equal')
plt.savefig('../../data/glamer/glamer_dim_scatter.png')
#plt.show()