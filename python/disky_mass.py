import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

radius=2.0 # kpc, radius of cylinder
length=60.0  # kpc, length of cylinder
res=100 # resolution on the circle

# setup square mass grid  
#mass=np.zeros((res+1,res+1,res*length/radius+1))

# setup distance grid

circle_x=np.linspace(0,radius,res+1)
circle_y=np.linspace(0,radius,res+1)

dx,dy=np.meshgrid(circle_x,circle_y)

# set up circle boolean array

circle_r=np.sqrt(dx**2+dy**2)

flag=np.less_equal(circle_r,radius)

cylin_d=np.linspace(0.001,length,res*length)
#cylin_d=np.logspace(-10,3,res*length)
#print cylin_d

#cylin_d=np.linspace(0.0001,60,100)

## define mass model


def disk(r,z):
	# use think here
	sig_d=209.5*10**6
	z_d=0.9
	Rd=3.31

#	if r<=Rd:
	rho_d=sig_d/2/z_d*np.exp(-z/z_d-r/Rd)
#	else:
#		rho_d=0.0

	return rho_d

def NFW(r):
	rho_h0=0.00846*10**9
	rh=20.2
	x=r/rh

	rho=rho_h0/x/(1+x)**2

	return rho

## calculate mass

Rd=3.31

disk_mass=0.0
halo_mass=0.0
dl=cylin_d[-1]-cylin_d[-2]


for i in range(len(cylin_d)):

	#dr=cylin_d[i]

	dr=np.sqrt(cylin_d[i]**2+circle_x**2) # R to center
	flag2=np.less_equal(dr,Rd) # flag for disk scale
	

	#if cylin_d[i]<=Rd:
	mass_gird=disk(dr,dy)*flag#*flag2
	mass=scipy.interpolate.RectBivariateSpline(circle_x,circle_y,mass_gird)
	mass=scipy.interpolate.RectBivariateSpline.integral(mass,0,radius,0,radius)
		#print mass

		#disk_mass=disk_mass+np.sum(mass_gird)*dl
	disk_mass=disk_mass+mass*dl
	#print 'disk_mass='
	#print disk_mass

	mass_gird=NFW(dr)*flag
	mass=scipy.interpolate.RectBivariateSpline(circle_x,circle_y,mass_gird)
	mass=scipy.interpolate.RectBivariateSpline.integral(mass,0,radius,0,radius)
	#halo_mass=halo_mass+np.sum(mass_gird)*dl
	halo_mass=halo_mass+mass*dl
	#print 'halo_mass='
	#print halo_mass

print disk_mass/(disk_mass+halo_mass)

'''

## plot profile

rho_disk=np.zeros(len(cylin_d))
rho_halo=np.zeros(len(cylin_d))

for i in range(len(cylin_d)):
	rho_disk[i]=disk(cylin_d[i],0)
	rho_halo[i]=NFW(cylin_d[i])


#print np.sum(rho_disk)/(np.sum(rho_disk)+np.sum(rho_halo))

plt.loglog(cylin_d,rho_disk,'r')
plt.loglog(cylin_d,rho_halo)

plt.show()

'''
