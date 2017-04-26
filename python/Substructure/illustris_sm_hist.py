import numpy as np
import matplotlib.pyplot as plt

## all lenses
tab = np.loadtxt('../../data/illustris_1/snap99_p1_all.txt')
all_ms, all_re = tab[:,0],tab[:,1]
tab = np.loadtxt('../../data/illustris_1/snap99_p2_all.txt')
all_ms, all_re= np.append(all_ms,tab[:,0]),np.append(all_re,tab[:,1])
tab = np.loadtxt('../../data/illustris_1/snap99_p3_all.txt')
all_ms, all_re= np.append(all_ms,tab[:,0]),np.append(all_re,tab[:,1])
# re mask
mask = all_re>=0.1
all_ms,all_re = all_ms[mask],all_re[mask]

## morphology
tab = np.loadtxt('../../data/illustris_1/snap99_p1_mor.txt')
mor_ms, mor_re = tab[:,0],tab[:,1]
tab = np.loadtxt('../../data/illustris_1/snap99_p2_mor.txt')
mor_ms, mor_re= np.append(mor_ms,tab[:,0]),np.append(mor_re,tab[:,1])
tab = np.loadtxt('../../data/illustris_1/snap99_p3_mor.txt')
mor_ms, mor_re= np.append(mor_ms,tab[:,0]),np.append(mor_re,tab[:,1])
# re mask
mask = mor_re>=0.1
mor_ms,mor_re = mor_ms[mask],mor_re[mask]

## kinematic
tab = np.loadtxt('../../data/illustris_1/snap99_p1_k.txt')
k_ms, k_re = tab[:,0],tab[:,1]
tab = np.loadtxt('../../data/illustris_1/snap99_p2_k.txt')
k_ms, k_re= np.append(k_ms,tab[:,0]),np.append(k_re,tab[:,1])
tab = np.loadtxt('../../data/illustris_1/snap99_p3_k.txt')
k_ms, k_re= np.append(k_ms,tab[:,0]),np.append(k_re,tab[:,1])
# re mask
mask = k_re>=0.1
k_ms,k_re = k_ms[mask],k_re[mask]

#bins = np.arange(9.8,11.8,0.2)
#all_his,b = np.histogram(np.log10(all_ms),bins)
#mor_his,b = np.histogram(np.log10(mor_ms),bins)
#k_his,b = np.histogram(np.log10(k_ms),bins)
bins=np.arange(0.,2.2,0.2)
all_his,b = np.histogram(all_re,bins)
mor_his,b = np.histogram(mor_re,bins)
k_his,b = np.histogram(k_re,bins)

print all_his,mor_his,k_his
mor_ra = mor_his.astype(float)/all_his.astype(float)
k_ra = k_his.astype(float)/all_his.astype(float)

#plt.hist(np.log10(k_ms),color = 'b',alpha=0.3,label='')
#plt.hist(np.log10(m_ms),color = 'b',label='morphology')
#plt.hist(np.log10(d_ms),color = 'r',label='kinematic')
plt.plot(bins[:-1]+0.2,mor_ra,'b',lw=2,linestyle='steps--',label='morphology')
plt.plot(bins[:-1]+0.2,k_ra,'r',lw=2,linestyle='steps',label='kinematic')
#plt.xlabel('log(M*)')
plt.xlabel('Einstein radius (")')
plt.ylabel('Disc system fraction')
#plt.xlim(10,11.5)
#plt.xlim(0,2.0)
plt.ylim(0,0.8)
plt.legend()
#plt.show()
plt.savefig('../../data/glamer/snap99_re_hist.png')