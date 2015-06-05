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

x=[0.0 , -0.0485 ,-1.1831,-1.7111]
y=[0.0,-0.1698,-0.6471 ,0.4626]

sx,sy=-6.471441e-01,1.721199e-01

fig=plt.figure()
plt.plot(x0,y0,'b+',ms=10,label='observe')
plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
plt.plot(sx,sy,'o',ms=5,mec='k',mfc='k')

plt.xlim(-1.8,0.2)
plt.ylim(-1.0,1.0)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.legend()
plt.show()

