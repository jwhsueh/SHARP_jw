import numpy as np
import matplotlib.pyplot as plt

tab = np.loadtxt('../../data/illustris_1/snap99_p1_kms.txt')
k_ms = tab
tab = np.loadtxt('../../data/illustris_1/snap99_p2_kms.txt')
k_ms = np.append(k_ms,tab)
tab = np.loadtxt('../../data/illustris_1/snap99_p3_kms.txt')
k_ms = np.append(k_ms,tab)

tab = np.loadtxt('../../data/illustris_1/snap99_p1_mms.txt')
m_ms = tab
tab = np.loadtxt('../../data/illustris_1/snap99_p2_mms.txt')
m_ms = np.append(m_ms,tab)
tab = np.loadtxt('../../data/illustris_1/snap99_p3_mms.txt')
m_ms = np.append(m_ms,tab)

tab = np.loadtxt('../../data/illustris_1/snap99_p1_dms.txt')
d_ms = tab
tab = np.loadtxt('../../data/illustris_1/snap99_p2_dms.txt')
d_ms = np.append(d_ms,tab)
tab = np.loadtxt('../../data/illustris_1/snap99_p3_dms.txt')
d_ms = np.append(d_ms,tab)

#plt.hist(np.log10(k_ms),color = 'b',alpha=0.3,label='')
plt.hist(np.log10(m_ms),color = 'b',label='morphology')
plt.hist(np.log10(d_ms),color = 'r',label='kinematic')
plt.xlabel('log(M*)')
plt.ylabel('Galaxy count')
plt.legend()
#plt.show()
plt.savefig('../../data/glamer/snap99_sm_hist.png')