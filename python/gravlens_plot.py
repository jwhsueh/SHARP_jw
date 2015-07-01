import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

## inputs for plotting

# gravlens file name

name='mode1.crit'

# lensed images
## system == B1555
# observe data
x0=[0.0, 0.0726,0.4117,0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

# model
x=[0.4108 , 0.1617 ,-0.00,0.0728]
y=[-0.0309,-0.3565,-0.0 ,0.0576]

# model paras
# source
sx,sy=2.011120e-01, -1.430824e-01

# mass profile centeroid
cx=[1.883924e-01,1.594156e-01] # x position
cy=[-1.607033e-01, -2.224752e-01] # y position


##----call functions in the end----##

def crit():
    ## read gravlens file
    t=np.loadtxt(name)

    ## critical curves
    x1,x2=t[:,0],t[:,4]
    y1,y2=t[:,1],t[:,5]

    n=x1.shape[0]
    Lx=np.zeros((n,2))
    Ly=np.zeros((n,2))

    for i in range(n):
        Lx[i,0],Lx[i,1]=x1[i],x2[i]
        Ly[i,0],Ly[i,1]=y1[i],y2[i]


    for i in range(n):
        plt.plot(Lx[i,:],Ly[i,:],'b')

#model_plot()

##----

def caus():
## caustics
    # read gravlens file
    t=np.loadtxt(name)

    u1,u2=t[:,2],t[:,6]
    v1,v2=t[:,3],t[:,7]

    n=u1.shape[0]
    Lu=np.zeros((n,2))
    Lv=np.zeros((n,2))

    for i in range(n):
        Lu[i,0],Lu[i,1]=u1[i],u2[i]
        Lv[i,0],Lv[i,1]=v1[i],v2[i]


    for i in range(n):
        plt.plot(Lu[i,:],Lv[i,:],'r--')

#model_plot()

##----

## souce & component positions
def model_plot():
    plt.plot(sx,sy,'o',ms=5,mec='k',mfc='k',label='source')
    plt.plot(cx,cy,'^',label='lens',mfc='k')

##----

def img_pos():

    plt.plot(x0,y0,'b+',ms=10,label='observe')
    plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
    
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')



##---- call functions from here ----##

crit()
caus()
model_plot()
img_pos()
plt.legend(loc=4)
plt.show()