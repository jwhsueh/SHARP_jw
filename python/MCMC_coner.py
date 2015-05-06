import numpy as np
import triangle
import matplotlib.pyplot as plt
#import matplotlib.patches as mpat

file_name=['SIE2_free_mcmc.dat','SIE2_free2_mcmc.dat','SIE2_free3_mcmc.dat','SIE2_free4_mcmc.dat']

k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp

    else:
        t=np.append(t,temp,axis=0)

    k=k+1



chi2=t[:,0]
sig1,sig2=t[:,1],t[:,6]
x1,x2=t[:,2],t[:,7]
y1,y2=t[:,3],t[:,8]
e1,e2=t[:,4],t[:,9]
pa1,pa2=t[:,5],t[:,10]

data=t[:,1:6]
print len(chi2)

# Generate some fake data.


#pick a parameter to plot
tot='total counts = %d' %(len(chi2))

## triangle plot

#ndim, nsamples = 10, len(chi2)

figure = triangle.corner(data, labels=[r"$\sigma$", r"$x$", r"$y$",
                                       r"$e$",r"$PA$"],
                         
                         quantiles=[0.16, 0.5, 0.84],
                         show_titles=True, title_args={"fontsize": 12})
figure.savefig("test.png")




