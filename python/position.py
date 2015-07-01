import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

## B1555
x0=[0.0, 0.0726,0.4117,0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

x=[0.4108 , 0.1617 ,-0.00,0.0728]
y=[-0.0309,-0.3565,-0.0 ,0.0576]

sx,sy=2.011120e-01, -1.430824e-01

## B0712
#x0=[0.0, -0.075,-1.185,-1.71]
#y0=[0.0, -0.16,-0.67,+0.46]

#x=[0.0 , -0.0583 ,-1.1751,-1.694]
#y=[0.0,-0.1693,-0.6757 ,0.4709]

#sx,sy=1.901373e-01, -1.433305e-01

#cx=[-8.494356e-01,-1.063329e+00]
#cy=[2.146172e-01, -2.157608e-03]

fig=plt.figure()
plt.plot(x0,y0,'b+',ms=10,label='observe')
plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
plt.plot(sx,sy,'o',ms=5,mec='k',mfc='k',label='source')
#plt.plot(cx,cy,'x',label='component')

#plt.xlim(-1.8,0.2)
#plt.ylim(-1.0,1.0)
plt.xlim(-0.1,0.5)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.legend()
plt.show()

