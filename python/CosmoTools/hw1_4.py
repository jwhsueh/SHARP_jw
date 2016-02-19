import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np

#parameters from WMAP 7yr result
global OM,OL,OK,h
OM=0.222+0.0449
OL=0.734
OK=1.0-OM-OL
h=0.71

#redshift array
z=np.linspace(0.00000001,10,10001)

#Friedmann E(z)
rEz=lambda x: (OM*(1.+x)**3+OK*(1+x)**2+OL)**(-0.5)

#comovinng distance
def DM(zi):
	dis=scipy.integrate.quad(rEz,0,zi)
	dis=3000.0/h*dis[0]
	return dis

##Luminosity distance
DL=np.zeros(10001)
for i in range(10001):
	DL[i]=(1.0+z[i])*DM(z[i])

#apparent magnitude
Mb=-20.6

am=np.zeros(10001)
am=Mb-5.0*(1.0-np.log10(DL)-6.0)

print am[0:10]

##Angular diameter distance
#DA=np.zeros(101)
#for i in range(101):
#	DA[i]=DM(z[i])/(1.0+z[i])

#fig1=plt.loglog(z,DL)
fig1=plt.semilogx(z,am)
plt.xlabel('redshift')
plt.ylabel('apparent magnitude')
plt.show()
