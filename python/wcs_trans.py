from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt

x = np.array([53,73,219,307])-256
y = np.array([224,323,480,223])-256

ps = 0.0058 #pix scale in arcsec

img = np.array([-1.0*x,y])*ps
print img[:,0]

theta = np.radians(130)
c,s = np.cos(theta),np.sin(theta)
R = np.array([[c,-1.0*s],[s,c]])
print R

A = np.zeros((2,4))
B = np.zeros((2,4))

for i in range(4):
	A[:,i] =  np.dot(R,img[:,i])

for i in range(4):
	B[:,i] = A[:,i]-A[:,1]

print B

#print img

plt.scatter(B[0],B[1])
plt.xlim(2.0,-1.0)
plt.ylim(-2.0,1.0)
plt.show()