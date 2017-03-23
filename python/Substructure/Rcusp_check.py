import numpy as np
import matplotlib.pyplot as plt


filepath='/Volumes/sting_1/data/Rcusp_c/'
#sie = np.loadtxt(filepath+'1798991_p3_64_Rcusp_sie.txt')
sie = np.loadtxt(filepath+'1798991_p3_64_Rcusp_ga.txt')
sie2 = np.loadtxt(filepath+'1798991_p3_128_Rcusp_ga.txt')
sie3 = np.loadtxt(filepath+'179899_p1_256_Rcusp_ga.txt')
pt = np.loadtxt(filepath+'179899_p1_64_Rcusp_s.txt')
ga = np.loadtxt(filepath+'179899_p1_64_Rcusp_ga.txt')
ga2= np.loadtxt(filepath+'179899_la_p1_64_Rcusp_ga.txt')
ga3= np.loadtxt(filepath+'179899_p1_256_Rcusp_ga.txt')


#sie2 = np.loadtxt(filepath+'177811_p2_64_Rcusp_sie.txt')
pt2 = np.loadtxt(filepath+'177811_p2_64_Rcusp_s.txt')
#ga2 = np.loadtxt(filepath+'177811_p2_64_Rcusp_ga.txt')
ga1282 = np.loadtxt(filepath+'177811_p2_128_Rcusp_ga.txt')

#sie3 = np.loadtxt(filepath+'177217_p3_64_Rcusp_sie.txt')
pt3 = np.loadtxt(filepath+'177217_p3_64_Rcusp_s.txt')
#ga3 = np.loadtxt(filepath+'177217_p3_64_Rcusp_ga.txt')
ga1283 = np.loadtxt(filepath+'177217_p3_128_Rcusp_ga.txt')
'''

sie_phi1 = sie[:,3]
#sie_phi1[sie_phi1>45] = 45-np.abs(45-sie_phi1[sie_phi1>45])
sie2_phi1 = sie2[:,3]
#sie2_phi1[sie2_phi1>45] = 45-np.abs(45-sie2_phi1[sie2_phi1>45])
sie3_phi1 = sie3[:,3]
#sie3_phi1[sie3_phi1>45] = 45-np.abs(45-sie3_phi1[sie3_phi1>45])

sie_phi0 = sie[:,2]
#sie_phi0[sie_phi0>120] = 120-np.abs(120-sie_phi0[sie_phi0>120])
'''
fig,ax=plt.subplots(2)
#plt.scatter(ga[:,2],np.abs(ga[:,1]),marker='o',color='b',label='MERLIN gauss fit, NN=64')
#plt.scatter(ga128[:,2],np.abs(ga128[:,1]),marker='*',color='r',label='NN=128')

#ax[0].scatter(sie[:,2],np.abs(sie[:,1]),marker='*',color='r',label='NN=64')
#ax[1].scatter(sie2[:,2],np.abs(sie2[:,1]),marker='*',color='b',label='NN=128')
#ax[2].scatter(sie3[:,2],np.abs(sie3[:,1]),marker='*',color='g',label='NN=256')

ax[0].scatter(ga[:,3],np.abs(ga[:,0]),marker='*',color='r',label='< 0.15*r')
ax[1].scatter(ga2[:,3],np.abs(ga2[:,0]),marker='*',color='b',label='< 0.25*r')
#ax[2].scatter(ga3[:,2],np.abs(ga3[:,1]),marker='*',color='g',label='NN=256')
#plt.scatter(pt[:,2],np.abs(pt[:,1]),marker='x',color='g',label='pt mag')

#plt.scatter(ga2[:,2],np.abs(ga2[:,1]),marker='o',color='b')
#plt.scatter(ga1282[:,2],np.abs(ga1282[:,1]),marker='*',color='r')
#plt.scatter(sie2[:,2],np.abs(sie2[:,1]),marker='*',color='r')
#plt.scatter(pt2[:,2],np.abs(pt2[:,1]),marker='x',color='g')

#plt.scatter(ga3[:,2],np.abs(ga3[:,1]),marker='o',color='b')
#plt.scatter(ga1283[:,2],np.abs(ga1283[:,1]),marker='*',color='r')
#plt.scatter(sie3[:,2],np.abs(sie3[:,1]),marker='*',color='r')
#plt.scatter(pt3[:,2],np.abs(pt3[:,1]),marker='x',color='g')
ax[0].set_title('Rcusp sanity check, NN=64')
plt.xlabel('delta phi')
plt.ylabel('|Rcusp|')
ax[0].legend(loc=2)
ax[1].legend(loc=2)
#ax[2].legend(loc=2)

#ax[0].set_xlim(40,140)
#ax[1].set_xlim(40,140)
#ax[2].set_xlim(40,140)
ax[0].set_ylim(0,0.3)
ax[1].set_ylim(0,0.3)
plt.ylim(0,0.3)
#plt.xlim(60,140)
plt.show()
#plt.savefig('../../data/glamer/Rcusp_sc_a.png')