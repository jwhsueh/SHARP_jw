import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

t=np.loadtxt('SIE2_mcmc_v3.dat')
n_para='e1'

chi2=t[:,0]
para1=t[:,1]

plt.figure(1)
n, bin, patches = plt.hist(para1,50)
plt.xlabel(n_para)
plt.ylabel('counts')

print len(chi2)

plt.figure(2)
plt.scatter(para1,chi2)
plt.xlabel(n_para)
plt.ylabel('chi^2')
#plt.legend()
plt.show()


