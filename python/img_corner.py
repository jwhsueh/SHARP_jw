import numpy as np
import triangle
import matplotlib.pyplot as plt
#import matplotlib.patches as mpat

file_name=['SIEsh_3_imgchain.dat']

k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp

    else:
        t=np.append(t,temp,axis=0)

    k=k+1


xa,ya,fa=t[:,0],t[:,1],t[:,2]
xb,yb,fb=t[:,3],t[:,4],t[:,5]
xc,yc,fc=t[:,6],t[:,7],t[:,8]
xd,yd,fd=t[:,9],t[:,10],t[:,11]

rab=fa/fb
rac=fa/fc
rad=fa/fd

data=t

## triangle plot

#ndim, nsamples = 10, len(chi2)

figure = triangle.corner(data, labels=[r"$x_A$", r"$y_A$",r"$f_A$", r"$x_B$",
                                       r"$y_B$",r"f_B",r"$x_C$",
                                       r"$y_C$",r"$f_C$" ,r"$x_D$", r"$y_D$",
                                    r"$f_D$"],
        #                         extents=[(120.,170.),(0.1,0.25),(-0.35,-0.2),(0.,0.5),(60.,110.),
                #                          (110.,160.),(0.1,0.25),(-0.3,-0.1),(0.75,0.95),(0.,20.)],

                         quantiles=[0.16, 0.5, 0.84],
                         show_titles=True, title_args={"fontsize": 12})
                         
figure.savefig("imgtest.png")




