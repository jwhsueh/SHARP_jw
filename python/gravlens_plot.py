import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

## inputs for plotting

# gravlens file name

name='B1555_expdisk_try5.crit'

## NOTE!! I've made change only on x axis of B1555

# lensed images

## system == B1555
# observe data
x0=[0.0, -0.0726,-0.4117,-0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

# model
x=[-0.4117 , -0.1629 ,-0.00,-0.0727]
y=[-0.0281,-0.3653,-0.00 ,0.0477]

# model paras
# source
sx,sy= -1.952576e-01, -1.493771e-01

# mass profile centeroid
cx=[ -1.818271e-01, -1.471605e-01] # x position
cy=[ -1.987580e-01, -2.056755e-01] # y position


'''
## system == B0712
# observe data
x0=[0.0, -0.075,-1.185,-1.71]
y0=[0.0, -0.16, -0.67,0.46]

# model
x=[0.0, -0.075,-1.185,-1.7103]
y=[0.0, -0.16, -0.67,0.4600]

# model paras
# source
sx,sy= -8.809485e-01,  2.514535e-01

# mass profile centeroid
cx=[ -8.335581e-01, -1.354868e+00] # x position
cy=[ 3.426335e-01, 1.572431e-01] # y position
'''


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
        Lx[i,0],Lx[i,1]=-x1[i],-x2[i]
        Ly[i,0],Ly[i,1]=y1[i],y2[i]

##


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
        Lu[i,0],Lu[i,1]=-u1[i],-u2[i]
        Lv[i,0],Lv[i,1]=v1[i],v2[i]

##


    for i in range(n):
        plt.plot(Lu[i,:],Lv[i,:],'r--')

#model_plot()

##----

## souce & component positions
def model_plot():
    plt.plot(sx,sy,'o',ms=5,mec='k',mfc='r',label='source')
    plt.plot(cx,cy,'^',label='lens',mfc='k')

##----

def img_pos():

    plt.plot(x0,y0,'b+',ms=10,label='observe')
    plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
    
    plt.xlabel(r'$\Delta  \alpha $ (arcsec)')
    plt.ylabel('$\Delta \delta$ (arcsec)')
    #plt.xlabel('arcsec')
    #plt.ylabel('arcsec')




##---- call functions from here ----##

crit()
caus()
model_plot()
img_pos()

plt.xlim(0.4,-0.8)
plt.ylim(-0.8,0.4)
plt.legend(loc=1)
plt.show()