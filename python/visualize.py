import numpy as np
import matplotlib.pyplot as plt

re=1024
L=1.5 #arcsec
pix=L/re #pixel size

x=np.linspace(-L/2.0,L/2.0,re)
y=np.linspace(-L/2.0,L/2.0,re)

xi,yi=np.meshgrid(x,y)

#parameters
b=0.2
e=0.56
q=1-e
esp=np.sqrt((1.0-q)/(1.0+q))

kappa=0.5*b/np.sqrt((1.0-esp)*xi**2+(1.0+esp)*yi**2)

fig1=plt.imshow(np.log(kappa))
#plt.clim(0.,50)
plt.colorbar()
plt.show()
