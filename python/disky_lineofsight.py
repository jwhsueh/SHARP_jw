import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import pyfits
import scipy.interpolate

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
#rho_h0=0.00846*10**9
rho_h0=0.00846*10**9
rh=20.2
dc=1#.6


# line of sight mass density
def disk(y,z,dtype):
	sigma=2*dtype.sig/2/dtype.zd*y*np.exp(-1*np.abs(z)/dtype.zd)*scipy.special.kn(1,y/dtype.Rd)

	return sigma

def NFW(y,z):
	x=np.sqrt(y**2+z**2)/rh
	
	#flag1=np.less_equal(x,1)
	#print flag1
	sigma1=2*rh*dc*rho_h0/(x**2-1)*(1-2/np.sqrt(1-x**2)*np.arctanh(np.sqrt((1-x)/(1+x))))
	sigma2=2*rh*dc*rho_h0/(x**2-1)*(1-2/np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1)/(1+x))))
	
	sigma1[x>1]=0.0

	sigma2[x<=1]=0.0

	sigma=sigma1+sigma2
	return sigma


def cylinder_massfraction(disk,halo,radius):
	# disk mass grid, halo mass grid, cylinder radius
	r=np.sqrt(y**2+z**2)
	flag=np.less_equal(r,radius)

	# flag mass
	disk=disk*flag
	halo=halo*flag

	# get rid of nan
	disk[:,0]=0.0
	halo[:,0]=0.0

	mass_d=scipy.interpolate.RectBivariateSpline(dz,dy,disk)
	mass_d=scipy.interpolate.RectBivariateSpline.integral(mass_d,-radius,radius,0,radius)

	mass_h=scipy.interpolate.RectBivariateSpline(dz,dy,halo)
	mass_h=scipy.interpolate.RectBivariateSpline.integral(mass_h,-radius,radius,0,radius)

	return mass_d/(mass_h+mass_d)



## plot line of sight desity profile
'''
# 1D plot
y=np.linspace(0,200,200*100)
y=np.array(y)
z=0.0
'''

## 2D plot
dy=np.linspace(0,10,10*100)
dz=np.linspace(-10,10,20*100)

y,z=np.meshgrid(dy,dz)

disk_thin=disk(y,z,thin)
disk_thick=disk(y,z,thick)
halo=NFW(y,z)

disk_mass=disk_thin+disk_thick
total_mass=disk_thin+disk_thick+halo
#total_mass=np.log10(total_mass)


'''
## cylinder flag

print cylinder_massfraction(disk_mass,halo,10)
print cylinder_massfraction(disk_mass,halo,4)
print cylinder_massfraction(disk_mass,halo,8)
print cylinder_massfraction(disk_mass,halo,10)
'''

hdu=pyfits.PrimaryHDU(disk_mass)
hdu.writeto('disk_project_mass.fits')

'''
## save figure

plt.imshow(total_mass,extent=[0,10,0,10])
plt.colorbar()
plt.xlabel('kpc')
plt.ylabel('kpc')
plt.title( r' $\Sigma_{total}$ (M_sun/kpc^2)')
#plt.show()
plt.savefig('total_project_mass.png',bbox_inches='tight')
'''

'''
## 1D projected mass profile

plt.loglog(y,halo,label='NFW')
plt.loglog(y,disk_thin,'r',label='disk_thin')
plt.loglog(y,disk_thick,'k',label='disk_thick')

plt.xlim(0,200)
plt.ylim(1e6,1e12)
plt.title('Edge-on view mass density profile, z=0')
plt.xlabel('r(kpc)')
plt.ylabel('Sigma')

plt.legend()

#plt.show()
plt.savefig('lineofsight_profile0.png',bbox_inches='tight')
'''

