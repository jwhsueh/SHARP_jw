import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter

fig, ((ax11,ax22),(ax1,ax2)) = plt.subplots(2,2,sharey=True,figsize=(6,6))

tab = np.loadtxt('../../models/snap99_179899/179899_area15_rcusp.txt')
tab2 = np.loadtxt('../../models/snap99_179899/179899_area25_rcusp.txt')
tab3 = np.loadtxt('../../models/snap99_179899/179899_area35_rcusp.txt')

rcusp,phi0 = tab[:,1],tab[:,2]
rcusp_2,phi0_2 = tab2[:,1],tab2[:,2]
rcusp_3,phi0_3 = tab3[:,1],tab3[:,2]

rfold,phi1 = tab[:,0],tab[:,3]
rfold_2,phi1_2 = tab2[:,0],tab2[:,3]
rfold_3,phi1_3 = tab3[:,0],tab3[:,3]

ax1.scatter(phi0_3,np.abs(rcusp_3),color='k',marker='o',label='e = 0.38')
ax1.scatter(phi0_2,np.abs(rcusp_2),color='r',marker='o')
ax1.scatter(phi0,np.abs(rcusp),color='b',marker='o')
ax1.set_ylim(0,0.7)
ax1.set_xlim(40,140)

ax11.scatter(phi1_3,np.abs(rfold_3),color='k',marker='o',label='e = 0.38')
ax11.scatter(phi1_2,np.abs(rfold_2),color='r',marker='o')
ax11.scatter(phi1,np.abs(rfold),color='b',marker='o')
ax11.set_ylim(0,0.7)
ax11.set_xlim(0,60)
#plt.legend(loc=2,scatterpoints=1)
#ax1.get_legend()

tab = np.loadtxt('../../models/snap99_179899/179899_d1_area15_rcusp.txt')
tab2 = np.loadtxt('../../models/snap99_179899/179899_d1_area25_rcusp.txt')
tab3 = np.loadtxt('../../models/snap99_179899/179899_d1_area35_rcusp.txt')

rcusp,phi0 = tab[:,1],tab[:,2]
rcusp_2,phi0_2 = tab2[:,1],tab2[:,2]
rcusp_3,phi0_3 = tab3[:,1],tab3[:,2]

rfold,phi1 = tab[:,0],tab[:,3]
rfold_2,phi1_2 = tab2[:,0],tab2[:,3]
rfold_3,phi1_3 = tab3[:,0],tab3[:,3]

ax2.scatter(phi0_3,np.abs(rcusp_3),color='k',marker='o',label='e = 0.38')
ax2.scatter(phi0_2,np.abs(rcusp_2),color='r',marker='o')
ax2.scatter(phi0,np.abs(rcusp),color='b',marker='o')
ax2.set_ylim(0,0.7)
ax2.set_xlim(40,140)

ax22.scatter(phi1_3,np.abs(rfold_3),color='k',marker='o',label='e = 0.28')
ax22.scatter(phi1_2,np.abs(rfold_2),color='r',marker='o')
ax22.scatter(phi1,np.abs(rfold),color='b',marker='o')
ax22.set_ylim(0,0.7)
ax22.set_xlim(0,60)
#plt.legend(loc=2,scatterpoints=1)
'''
tab = np.loadtxt('../../models/snap99_179899/179899_d3_area15_rcusp.txt')
tab2 = np.loadtxt('../../models/snap99_179899/179899_d3_area25_rcusp.txt')
tab3 = np.loadtxt('../../models/snap99_179899/179899_d3_area35_rcusp.txt')

rcusp,phi0 = tab[:,1],tab[:,2]
rcusp_2,phi0_2 = tab2[:,1],tab2[:,2]
rcusp_3,phi0_3 = tab3[:,1],tab3[:,2]

rfold,phi1 = tab[:,0],tab[:,3]
rfold_2,phi1_2 = tab2[:,0],tab2[:,3]
rfold_3,phi1_3 = tab3[:,0],tab3[:,3]

ax33.scatter(phi0_3,np.abs(rcusp_3),color='k',marker='o',label='e = 0.38')
ax33.scatter(phi0_2,np.abs(rcusp_2),color='r',marker='o')
ax33.scatter(phi0,np.abs(rcusp),color='b',marker='o')
ax33.set_ylim(0,0.7)
ax33.set_xlim(40,140)

ax3.scatter(phi1_3,np.abs(rfold_3),color='k',marker='o',label='e = 0.18')
ax3.scatter(phi1_2,np.abs(rfold_2),color='r',marker='o')
ax3.scatter(phi1,np.abs(rfold),color='b',marker='o')
ax3.set_ylim(0,0.7)
ax3.set_xlim(0,60)
#plt.legend(loc=2,scatterpoints=1)
'''
#fig.subplots_adjust(wspace=0)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],visible=False)
fig.text(0.06,0.25,'|Rcusp|',va='center',rotation='vertical')
fig.text(0.06,0.75,'|Rfold|',va='center',rotation='vertical')
fig.text(0.5,0.03,'delta phi',ha='center')
fig.text(0.5,0.48,'phi 1',ha='center')

fig.text(0.16,0.85,'SIE',ha='left')
fig.text(0.16,0.82,'e=0.38',ha='left')
fig.text(0.58,0.85,'SIE+ExpDisc',ha='left')
fig.text(0.58,0.82,'e=0.76',ha='left')
fig.text(0.58,0.79,'Disc Mass=15%',ha='left')
#fig.text(0.72,0.85,'e=0.18',ha='center')

'''
H,xbin,ybin = np.histogram2d(rfold,phi1,bins=(np.linspace(0,0.5,30),np.linspace(40,140,30)))
H = gaussian_filter(H,1.0)
H2,xbin,ybin = np.histogram2d(rfold_2,phi1_2,bins=(np.linspace(0,0.5,30),np.linspace(40,140,30)))
H2 = gaussian_filter(H2,1.0)
#plt.imshow(H2,origin='lower',cmap=cm.Blues,extent=[ybin2[0],ybin2[-1],xbin2[0],xbin2[-1]])
#plt.imshow(H,origin='lower',cmap=cm.Greys,extent=[ybin2[0],ybin2[-1],xbin2[0],xbin2[-1]])
#plt.gca().set_aspect(ybin2[-1]/xbin2[-1])
#plt.colorbar()
#plt.ylim(0.05,0.5)

plt.contour(H2,colors='k',extent=[ybin[0],ybin[-1],xbin[0],xbin[-1]])
plt.contour(H,colors='r',extent=[ybin[0],ybin[-1],xbin[0],xbin[-1]])


plt.scatter(phi1_3,rfold_3,color='b',marker='o',label='r < 0.35*rc')
plt.scatter(phi1_2,rfold_2,color='r',marker='o',label='r < 0.25*rc')
plt.scatter(phi1,rfold,color='k',marker='o',label='r < 0.15*rc')

plt.title('SIE, e = 0.38')
plt.xlabel('phi 1')
plt.ylabel('|Rfold|')
plt.ylim(0,0.3)
plt.xlim(0,70)
#plt.xlim(20,180)
plt.legend(loc=2,scatterpoints=1)
'''
## subplots



plt.show()
#plt.savefig('../../data/gravlens_pt/elp_e38_rfold.png')