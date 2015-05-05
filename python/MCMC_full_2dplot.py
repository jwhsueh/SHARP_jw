import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

t=np.loadtxt('SIE2_free_mcmc.dat')


chi2=t[:,0]
sig1,sig2=t[:,1],t[:,6]
x1,x2=t[:,2],t[:,7]
y1,y2=t[:,3],t[:,8]
e1,e2=t[:,4],t[:,9]
pa1,pa2=t[:,5],t[:,10]

#pick 2 parameters to plot
para1=sig1
para2=e1
n_para='sigma1'
tot='total counts = %d' %(len(chi2))

plt.figure(1)
n, bin, patches = plt.hist(para1,50,label=tot)
plt.xlabel(n_para)
plt.ylabel('counts')

print len(chi2)

#plt.figure(2)
#plt.scatter(para1,chi2)
#plt.xlabel(n_para)
#plt.ylabel('chi^2')
plt.legend()
plt.show()


