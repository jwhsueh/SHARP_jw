import numpy as np
import matplotlib.pyplot as plt

data_path = '/Users/jwhsueh/Documents/SHARP_jw/data/glamer/'

edge = np.loadtxt(data_path+'snap99_sie_cusp_bin.txt')
#print edge

#edge = np.loadtxt(data_path+'snap99_sie_fold_bin.txt')#+2.5
#print edge

sie = np.loadtxt(data_path+'snap99_sie_cusp_pd.txt') 
sie_sig = np.loadtxt(data_path+'snap99_sie_cusp_sig.txt') 
#sie = np.transpose(sie)
print sie.shape

com = np.loadtxt(data_path+'snap99_elp_cusp_pd.txt') 
com_sig = np.loadtxt(data_path+'snap99_elp_cusp_sig.txt') 

edg = np.loadtxt(data_path+'snap99_edg_fold_pd.txt') 

fa = np.loadtxt(data_path+'snap99_fa_fold_pd.txt') 
#com = np.transpose(com)

plt.plot(edge[2:],sie[4,2:],color='k',linestyle='-.',label='SIE',lw=2)
plt.plot(edge[:],sie[3,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[2,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[1,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[0,:],color='k',linestyle='-.',lw=2)
plt.fill_between(edge[2:],sie[4,2:]-sie_sig[4,2:],sie[4,2:]+sie_sig[4,2:],color='k',alpha=0.1)
plt.fill_between(edge[:],sie[2,:]-sie_sig[2,:],sie[2,:]+sie_sig[2,:],color='k',alpha=0.1)
plt.fill_between(edge[0:],sie[3,:]-sie_sig[3,:],sie[3,:]+sie_sig[3,:],color='k',alpha=0.3)
plt.fill_between(edge[0:],sie[0,:]-sie_sig[0,:],sie[0,:]+sie_sig[0,:],color='k',alpha=0.1)
plt.fill_between(edge[0:],sie[1,:]-sie_sig[1,:],sie[1,:]+sie_sig[1,:],color='k',alpha=0.3)

plt.plot(edge[:],com[4,:],color='r',linestyle='-.',label='ell',lw=2)
plt.plot(edge,com[3,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[2,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[1,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[0,:],color='r',linestyle='-.',lw=2)
plt.fill_between(edge[:],com[4,:]-com_sig[4,:],com[4,:]+com_sig[4,:],color='r',alpha=0.1)
plt.fill_between(edge[0:],com[2,:]-com_sig[2,:],com[2,:]+com_sig[2,:],color='r',alpha=0.1)
plt.fill_between(edge[0:],com[3,:]-com_sig[3,:],com[3,:]+com_sig[3,:],color='r',alpha=0.3)
plt.fill_between(edge[0:],com[0,:]-com_sig[0,:],com[0,:]+com_sig[0,:],color='r',alpha=0.1)
plt.fill_between(edge[0:],com[1,:]-com_sig[1,:],com[1,:]+com_sig[1,:],color='r',alpha=0.3)
'''
plt.plot(edge[:],edg[4,:],color='b',linestyle='-',label='edge-on',lw=2)
plt.plot(edge,edg[3,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[2,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[1,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[0,:],color='b',linestyle='-',lw=2)
'''
'''
plt.plot(edge[:],fa[4,:],color='g',linestyle='-',label='face-on')
#plt.plot(edge,fa[3,:],color='g',linestyle='-')
plt.plot(edge,fa[2,:],color='g',linestyle='-')
#plt.plot(edge,fa[1,:],color='g',linestyle='-')
plt.plot(edge,fa[0,:],color='g',linestyle='-')
'''
#print com[0,:]

plt.legend(loc=2)
plt.xlim(40,110)
#plt.xlim(5,40)
plt.ylim(0,0.7)

#plt.xlabel('phi 1')
#plt.ylabel('|R_fold|')

plt.xlabel('delta phi')
plt.ylabel('|R_cusp|')
#plt.title('1%, 10%, 50% contour')

#plt.show()
plt.savefig(data_path+'snap99_elp_com_c.png')

