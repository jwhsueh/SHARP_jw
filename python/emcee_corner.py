import numpy as np
import triangle
import matplotlib.pyplot as plt
#import matplotlib.patches as mpat

file_name=['exp_chain.dat']

k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp

    else:
        t=np.append(t,temp,axis=0)

    k=k+1



#chi2=t[:,0]
#sig1,sig2=t[:,1],t[:,6]
#x1,x2=t[:,2],t[:,7]
#y1,y2=t[:,3],t[:,8]
#e1,e2=t[:,4],t[:,9]
#pa1,pa2=t[:,5],t[:,10]

#data=t[:,:11] #ignore source posotions first
#data=t[:,:9]
data=t
#print len(chi2)


## triangle plot

#ndim, nsamples = 10, len(chi2)

figure = triangle.corner(data, labels=[r"$b1$", r"$x1$", r"$y1$",
                                       r"$e1$",r"$PA1$",
                                       r"$b2$", r"$x2$", r"$y2$",
                                       r"$e2$",r"$PA2$",r"$rs$"
                                       ,r"x_s",r"y_s"],
#                         extents=[(120.,170.),(0.1,0.25),(-0.35,-0.2),(0.,0.5),(60.,110.),
#                                  (110.,160.),(0.1,0.25),(-0.3,-0.1),(0.75,0.95),(0.,20.)],
                         #range=[0.999,0.999,0.999,0.999,0.999,
                         #       0.999,0.999,0.999,0.999,0.999],
                         quantiles=[0.16, 0.5, 0.84],
                         show_titles=True, title_args={"fontsize": 12})
                         
figure.savefig("test.png")




