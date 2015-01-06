import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

x0=[0.0, -0.0726,-0.4117,-0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

x=[-0.4120, -0.1580 ,0.0017,-0.0734]
y=[-0.0278,-0.3608  ,0.0022 ,0.0456]

sx,sy=-2.036227e-01,-1.644783e-01

fig=plt.figure()
plt.plot(x0,y0,'b+',ms=10,label='observe')
plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='model')
plt.plot(sx,sy,'o',ms=5,mec='k',mfc='k')

plt.xlim(0.2,-0.6)
plt.ylim(-0.5,0.2)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
#plt.legend()
plt.show()

