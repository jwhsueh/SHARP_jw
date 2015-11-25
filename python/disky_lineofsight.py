import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import pyfits
import scipy.interpolate
import scipy.integrate

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

# bulge
alpha=1.8
r0=0.075
rcut=2.1
q=0.5
rho_b0=9.93e10


# line of sight mass surface density
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

## bulge (3D profile)

def cylinder_massfraction(disk,bulge,halo,radius):
	# disk mass grid, bulge mass grid, halo mass grid, cylinder radius
	# output: log(M_disk),log(M_bulge),log(M_halo),log(M_total),disk mass fraction

	r=np.sqrt(y**2+z**2)
	flag=np.less_equal(r,radius)

	# flag mass
	disk=disk*flag
	halo=halo*flag
	bulge=bulge*flag

	# get rid of nan
	disk[:,0]=0.0
	halo[:,0]=0.0

	mass_d=scipy.interpolate.RectBivariateSpline(dz,dy,disk)
	mass_d=scipy.interpolate.RectBivariateSpline.integral(mass_d,-radius,radius,0,radius)

	mass_b=scipy.interpolate.RectBivariateSpline(dz,dy,bulge)
	mass_b=scipy.interpolate.RectBivariateSpline.integral(mass_b,-radius,radius,0,radius)

	mass_h=scipy.interpolate.RectBivariateSpline(dz,dy,halo)
	mass_h=scipy.interpolate.RectBivariateSpline.integral(mass_h,-radius,radius,0,radius)

	return np.log10(mass_d),np.log10(mass_b), np.log10(mass_h), np.log10(mass_h+mass_d+mass_b),mass_d/(mass_h+mass_d+mass_b)
	#return mass_d/(mass_h+mass_d)

def num_bulge(y,z):

	sigma=np.zeros((len(z),len(y)))

	for i in range(len(z)): # vertical index
		print i
		for j in range(len(y)):

			bulge3d=lambda R: R*np.exp(-(R**2+z[i]**2/q**2))/(np.sqrt(R**2-y[j]**2)*(1+np.sqrt(R**2+z[i]**2/q**2)/r0)**alpha)
			I=scipy.integrate.quad(bulge3d,y[j],np.inf)

			sigma[i,j]=2*rho_b0*I[0]

	return sigma


## plot line of sight desity profile
'''
# 1D plot
y=np.linspace(0,200,200*100)
y=np.array(y)
z=0.0
'''

## 2D plot
radius=30 # kpc
dy=np.linspace(0,radius,radius*10)
dz=np.linspace(-1*radius,radius,2*radius*10)

y,z=np.meshgrid(dy,dz)

disk_thin=disk(y,z,thin)
disk_thick=disk(y,z,thick)
halo=NFW(y,z)
bulge=num_bulge(dy,dz) # numerical bulge doesn't take meshgrid

disk_mass=disk_thin+disk_thick
total_mass=disk_thin+disk_thick+halo+bulge
#total_mass=np.log10(total_mass)



## cylinder flag
'''
enclose_mass=np.zeros((5,4))
enclose_mass[:,0]=cylinder_massfraction(disk_mass,halo,2)
enclose_mass[:,1]=cylinder_massfraction(disk_mass,halo,4)
enclose_mass[:,2]=cylinder_massfraction(disk_mass,halo,8)
enclose_mass[:,3]=cylinder_massfraction(disk_mass,halo,10)

plt.scatter(np.array([2,4,8,10]),enclose_mass[2,:],label='disk+halo',marker='^')
plt.scatter(np.array([2,4,8,10]),enclose_mass[0,:],label='disk',color='r')
plt.scatter(np.array([2,4,8,10]),enclose_mass[1,:],color='g',label='halo')
plt.xlabel('kpc')
plt.ylabel('log(M)')
plt.title('log(enclosed total mass)')
plt.legend(loc=4)

'''
mass_r=np.array([1,2,4,6,8,10,15,20,25,30])
mass_f=np.zeros((len(mass_r),5))

for i in range(len(mass_r)):
	mass_f[i,:]=cylinder_massfraction(disk_mass,bulge,halo,mass_r[i])

#print mass_f.shape

plt.figure(1)
plt.scatter(mass_r,mass_f[:,4])
plt.axvline(x=2,ymax=0.7)
plt.xlabel('kpc')
plt.ylabel('mass fraction')
plt.title('disk mass fraction')

#plt.show()	
plt.savefig('disk mass fraction.png',bbox_inches='tight')


plt.figure(2)
plt.scatter(mass_r,mass_f[:,3],label='total',marker='^')
plt.scatter(mass_r,mass_f[:,0],label='disk',color='r')
plt.scatter(mass_r,mass_f[:,2],color='g',label='halo')
plt.scatter(mass_r,mass_f[:,1],label='bulge')
plt.axvline(x=2,ymax=11.5)
plt.xlabel('kpc')
plt.ylabel('log(M)')
plt.title('log(enclosed total mass)')
plt.legend(loc=2)

plt.savefig('total enclose mass.png',bbox_inches='tight')

#hdu=pyfits.PrimaryHDU(disk_mass)
#hdu.writeto('disk_project_mass.fits')


## save figure
plt.figure(3)
plt.imshow(total_mass,extent=[0,radius,0,radius])
plt.colorbar()
plt.xlabel('kpc')
plt.ylabel('kpc')
plt.title( r' $\Sigma_{total}$ (M_sun/kpc^2)')
#plt.show()
plt.savefig('total_project_mass.png',bbox_inches='tight')


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

