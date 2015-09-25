import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

## inputs for plotting

# input files

name='/Users/jwhsueh/Documents/SHARP_jw/DK02/mock.crit'

name_findimg='/Users/jwhsueh/Documents/SHARP_jw/DK02/gravlens_mock_findimg.dat'

# glafic

#name_opt='exp_optresult.dat'

# number of components
ncom=2

## NOTE!! I've made change only on x axis of B1555

global x,y,x0,y0,sx,sy,cx,cy

# lensed images
'''
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
cx=[ -8.7e-01, -1.08] # x position
cy=[ 0.15, 0.08] # y position
'''

##----read model details from files [glafic]----##

def read_lens():
	## read glafic lens optimal result
	t=open(name_opt,'r')

        optfile=[]

        for line in t.readlines():
            optfile.append(line)

        src=optfile[-2]
        comp1=optfile[-4]
        comp2=optfile[-3]

        src=src.split()
        comp1=comp1.split()
        comp2=comp2.split()
    
        sx=float(src[2])
        sy=float(src[3])
        cx=[float(comp1[3]),float(comp2[3])]
        cy=[float(comp1[4]),float(comp2[4])]

	return sx,sy,cx,cy

##----read findimg positions from files [gravlens]----##

def read_findimg():
	## read gravlens image positions
	t=np.loadtxt(name_findimg)

    	x=t[:,0]
    	y=t[:,1]

	return x,y


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
        #        Lx[i,0],Lx[i,1]=-x1[i],-x2[i]
        Lx[i,0],Lx[i,1]=x1[i],x2[i]
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
        #        Lu[i,0],Lu[i,1]=-u1[i],-u2[i]
        Lu[i,0],Lu[i,1]=u1[i],u2[i]
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

#sx,sy,cx,cy=read_lens()
x,y=read_findimg()

plt.figure(figsize=(5.7,5.7))

crit()
caus()
#model_plot()
x0,y0=x,y
print x,y
img_pos()

#plt.xlim(0.4,-0.8)
#plt.ylim(-0.8,0.4)
#plt.legend(loc=1)

plt.show()


