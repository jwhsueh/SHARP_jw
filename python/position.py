import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

## B1555
#x0=[0.0, -0.0726,-0.4117,-0.1619]
#y0=[0.0, 0.0480, -0.0280,-0.3680]

#x=[-0.4117 , -0.1743 ,0.0005,-0.0727]
#y=[-0.0283,-0.3628,0.0005 ,0.0477]

#sx,sy=-1.890538e-01,-1.449221e-01

## B0712
x0=[0.0, -0.075,-1.185,-1.71]
y0=[0.0, -0.16,-0.67,+0.46]

x=[0.0 , -0.0466 ,-1.1768,-1.7049]
y=[0.0,-0.1747,-0.6860 ,0.4673]

sx,sy=-8.754612e-01,  1.457479e-01

cx=[-8.865909e-01,-8.558520e-01]
cy=[ 9.303342e-02,  4.755138e-01]

fig=plt.figure()
plt.plot(x0,y0,'b+',ms=10,label='observe')
plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
plt.plot(sx,sy,'o',ms=5,mec='k',mfc='k',label='source')
#plt.plot(cx,cy,'x',label='component')

plt.xlim(-1.8,0.2)
plt.ylim(-1.0,1.0)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.legend()
plt.show()

