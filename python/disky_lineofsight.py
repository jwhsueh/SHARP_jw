import numpy as np
import matplotlib.pyplot as plt
import scipy.special

# thin disk
class thin:
	zd=0.3
	Rd=2.9
	sig=816.6*10**6

# thick disk
class thick:
	zd=0.9
	Rd=3.31
	sig=209.5*10**6

# NFW
rho_h0=0.00846*10**9
rh=20.2
dc=140  # current critical density, M_sun/kpc^3





# line of sight mass density
def disk(y,z,dtype):
	sigma=2*dtype.sig/2/dtype.zd*y*np.exp(-z/dtype.zd)*scipy.special.kn(1,y/dtype.Rd)

	return sigma

def NFW(y,z):
	x=np.sqrt(y**2+z**2)/rh
	
	flag1=np.less_equal(x,1)
	sigma1=2*rh*dc*rho_h0/(x**2-1)*(1-2/np.sqrt(1-x**2)*np.arctanh(np.sqrt((1-x)/(1+x))))*flag1
	#elif x==1:
	#	sigma=2*rh*dc*rho_h0/3
	flag2=np.greater(x,1)
	sigma2=2*rh*dc*rho_h0/(x**2-1)*(1-2/np.sqrt(1-x**2)*np.arctan(np.sqrt((1-x)/(1+x))))*flag2

	sigma=sigma1+sigma2
	return sigma

## plot line of sight desity profile

y=np.linspace(0,200,200*100)
y=np.array(y)
z=0.0

halo=np.zeros(len(y))

disk_thin=disk(y,z,thin)
disk_thick=disk(y,z,thick)
halo=NFW(y,z)
#print y
#print halo

plt.loglog(y,halo,label='NFW')
plt.loglog(y,disk_thin,'r',label='disk_thin')
plt.loglog(y,disk_thick,'k',label='disk_thick')

plt.xlim(0,200)
plt.ylim(1e6,1e12)
plt.title('Edge-on view mass density profile, z=0')
plt.xlabel('r(kpc)')
plt.ylabel('Mass(M_sun)')

plt.legend()

#plt.show()
plt.savefig('lineofsight_profile0.png',bbox_inches='tight')

