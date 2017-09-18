import numpy as np
import matplotlib.pyplot as plt

path = '../../data/sub_gravlens/'
table = np.loadtxt(path+'/B1422_realization1/flux_result.txt')
#table = table[:1500,:]
fluxA,fluxB,fluxC,fluxD = table[:,2],table[:,1],table[:,0],table[:,3]
chi2 = table[:,5]
mod_fa = np.array([fluxB/fluxA,fluxC/fluxA,fluxD/fluxA])
print len(mod_fa[0])

fsub = np.loadtxt(path+'/B1422_realization1/real_fsub.txt')
mod_fsub = fsub[:5000]

## --- percentage threshold
p = np.linspace(0.01,0.2,20)
idx = 5000.*p
print idx
std_array = np.zeros(len(p))
chi2_sort = np.sort(chi2)

## --- absolute threshold
th = np.linspace(1,20,20)

for i in range(len(p)):
	th = chi2_sort[idx[i]]
	print th
	mask = chi2<th

	std_array[i] = np.std(mod_fsub[mask])
'''
plt.plot(p,std_array)
#plt.plot(th,std_array)
plt.title('fsub standard deviation')
#plt.xlabel('$\epsilon$ (percent)')
plt.xlabel('$\epsilon$ (absolute)')
#plt.show()
#plt.savefig('mock_abs_th.png')
'''
th = chi2_sort[idx[0]]
#th = 8.0
mask = chi2<1.0
print chi2[mask]
print mod_fsub[mask]

mock_fa = np.array([0.83748,0.47,0.03])
#mock_fsub = np.zeros(lens)
mock_fsub = 0.01

plt.scatter(mod_fsub,mod_fa[0],marker='.',alpha=0.5,label='random sample')
#plt.scatter(mod_fsub[mask],mod_fa[2][mask],marker='o',color='b',label='best 9%')

plt.scatter(mod_fsub[mask],mod_fa[0][mask],marker='o',color='r',label='chi2<8')
plt.scatter(mock_fsub,mock_fa[0],marker='*',facecolor='r',edgecolor='k',s=150,label='true')
#print len(mod_fsub),len(np.exp(-1*chi2**2)*5.0)
#plt.hist(mod_fsub,bins=np.logspace(-3,-1,20),weights=np.exp(-1*chi2/2))

plt.xscale('log')
#plt.plot([0.001,0.05],[mock_fa[2]*0.9,mock_fa[2]*0.9],linestyle='--',color='r',alpha=0.7,lw=2)
#plt.plot([0.001,0.05],[mock_fa[2]*1.1,mock_fa[2]*1.1],linestyle='--',color='r',alpha=0.7,lw=2)
#plt.plot(np.array([0.001,0.05]),np.array([0.83*1.1,0.83*1.1]))
plt.xlim(0.001,0.05)
#plt.ylim(mock_fa[0]+0.018,mock_fa[0]-0.018)
plt.xlabel('fsub')
#plt.ylabel('fD/fA')
#plt.title('flux ratio D/A')
#plt.legend(loc=3,scatterpoints=1)
plt.show()
#plt.savefig('mock_DA.png')
