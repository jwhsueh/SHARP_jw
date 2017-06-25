import numpy as np
import matplotlib.pyplot as plt

data_path = '/Users/jwhsueh/Documents/SHARP_jw/data/glamer/'

#edge = np.loadtxt(data_path+'snap99_sie_cusp_bin.txt')
#print edge

edge = np.loadtxt(data_path+'snap99_sie_fold_bin.txt')-2.5
print edge

sie = np.loadtxt(data_path+'snap99_sie_fold_pd.txt') 
#sie = np.transpose(sie)
print sie.shape

com = np.loadtxt(data_path+'snap99_elp_fold_pd.txt') 
#com = np.transpose(com)

plt.plot(edge[:],sie[4,:],color='k',linestyle='-.',label='SIE')
plt.plot(edge,sie[3,:],color='k',linestyle='-.')
plt.plot(edge,sie[2,:],color='k',linestyle='-.')
plt.plot(edge,sie[1,:],color='k',linestyle='-.')
plt.plot(edge,sie[0,:],color='k',linestyle='-.')

plt.plot(edge[:],com[4,:],color='r',linestyle='-',label='ell')
plt.plot(edge,com[3,:],color='r',linestyle='-')
plt.plot(edge,com[2,:],color='r',linestyle='-')
plt.plot(edge,com[1,:],color='r',linestyle='-')
plt.plot(edge,com[0,:],color='r',linestyle='-')

plt.legend(loc=2)
#plt.xlim(40,110)
plt.xlim(5,40)
plt.ylim(0,0.7)

plt.xlabel('Delta phi')
plt.ylabel('|R_fold|')

#plt.show()
plt.savefig(data_path+'snap99_elp_com_f.png')

