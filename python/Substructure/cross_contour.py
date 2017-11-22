import numpy as np
import matplotlib.pyplot as plt

data_path = '/Users/jwhsueh/Documents/SHARP_jw/data/glamer/'

edge = np.loadtxt(data_path+'snap99_sie_fold_bin.txt')
#print edge

#edge = np.loadtxt(data_path+'snap99_sie_fold_bin.txt')#+2.5
#print edge

sie = np.loadtxt(data_path+'snap99_sie_fold_pd.txt') 
#sie = np.transpose(sie)
print sie.shape

com = np.loadtxt(data_path+'snap99_elp_fold_pd.txt') 

edg = np.loadtxt(data_path+'snap99_edg_fold_pd.txt') 
fa = np.loadtxt(data_path+'snap99_fa_fold_pd.txt') 
#com = np.transpose(com)
'''
plt.plot(edge[2:],sie[4,2:],color='k',linestyle='-.',label='SIE',lw=2)
plt.plot(edge,sie[3,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[2,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[1,:],color='k',linestyle='-.',lw=2)
plt.plot(edge,sie[0,:],color='k',linestyle='-.',lw=2)
'''
plt.plot(edge[:],com[4,:],color='r',linestyle='-.',label='ell',lw=2)
plt.plot(edge,com[3,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[2,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[1,:],color='r',linestyle='-.',lw=2)
plt.plot(edge,com[0,:],color='r',linestyle='-.',lw=2)

plt.plot(edge[:],edg[4,:],color='b',linestyle='-',label='edge-on',lw=2)
plt.plot(edge,edg[3,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[2,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[1,:],color='b',linestyle='-',lw=2)
plt.plot(edge,edg[0,:],color='b',linestyle='-',lw=2)
'''
plt.plot(edge[:],fa[4,:],color='g',linestyle='-',label='face-on')
#plt.plot(edge,fa[3,:],color='g',linestyle='-')
plt.plot(edge,fa[2,:],color='g',linestyle='-')
#plt.plot(edge,fa[1,:],color='g',linestyle='-')
plt.plot(edge,fa[0,:],color='g',linestyle='-')
'''
#print com[0,:]

plt.legend(loc=2)
#plt.xlim(40,110)
plt.xlim(5,40)
plt.ylim(0,0.7)

plt.xlabel('phi 1')
plt.ylabel('|R_fold|')

#plt.xlabel('delta phi')
#plt.ylabel('|R_cusp|')
#plt.title('1%, 10%, 50% contour')

#plt.show()
plt.savefig(data_path+'snap99_ee_com_f.png')

