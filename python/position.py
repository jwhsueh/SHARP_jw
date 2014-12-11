import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat

x0=[0.0, 0.0726,0.4117,0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

#SIE_best
#x=[-0.0194,0.0595,0.4332,0.1378]
#y=[0.0654,0.0979,0.0150,-0.6783]

#x=[4.279748e-01,1.493585e-03,7.116276e-02,1.685729e-01]
#y=[-3.424671e-02,1.355427e-02 ,5.382633e-02,-3.834915e-01 ]

x=[0.4116,0.1638,0.0724,0.0000]
y=[-0.0283,-0.3651  ,0.0489 ,0.0000 ]

fig=plt.figure()
plt.plot(x0,y0,'b.',label='observe')
plt.plot(x,y,'r.',label='model')

plt.xlim(-0.6,0.2)
plt.ylim(-0.6,0.2)
plt.legend()
plt.show()

