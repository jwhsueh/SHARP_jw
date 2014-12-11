import numpy as np
import matplotlib.pyplot as plt

x0=[0.0, 0.0726,0.4117,0.1619]
y0=[0.0, 0.0480, -0.0280,-0.3680]

x=[]
y=[]

with open('./gravlens/extend1.txt','r') as p:
   i=0
   for line in p:
        if i>1:
	    s=[]
	    s.append(line.split())
	    num=s[0]
            x.append(float(num[0]))
            y.append(float(num[1]))
        else:
            i=i+1

x=np.array(x)
y=np.array(y)

x=-1.0*x

fig=plt.figure()
plt.plot(x0,y0,'b.',label='observe')
plt.plot(x,y,'r.',label='model')
plt.legend()
plt.show()
