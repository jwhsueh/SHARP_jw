import numpy as np
import triangle
import matplotlib.pyplot as plt
#import matplotlib.patches as mpat

file_name=['SIE2_free_mcmc.dat','SIE2_free2_mcmc.dat','SIE2_free3_mcmc.dat','SIE2_free4_mcmc.dat',
           'SIE2_free5_mcmc.dat','SIE2_free6_mcmc.dat']

# burn-out stage setting
b=50

k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp[b:,:]
    
    else:
        t=np.append(t,temp[b:,:],axis=0)
    
    k=k+1



chi2=t[:,0]
data=t[:,1:] # parameters
ndim=10

s=np.zeros(ndim)
for i in range(ndim):
    s[i]=np.median(data[:,i])
    #print np.median(data[:,i])

print s