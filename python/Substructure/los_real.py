import numpy as np
import scipy.stats as st
from astropy import cosmology
from scipy.integrate import trapz
import NFWprofile as NFW
import sys

### ---- calculate the total mass within z-range
lens='B1422'
zl = 0.34
zs=3.62
Re_arc = 0.75 # arcsec
Re = np.radians(Re_arc/3600.)*cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h #Mpc
#print Re

f_low,f_hi = 0.8,2.0

m_c = 0.0

path = '/Volumes/sting_1/subs/'+lens
real_folder = '/los_00/'

n_real =8000
n_start=int(sys.argv[1])

def massf(z,mass):
	a = 0.7689
	p = 0.2536
	AA = 0.3295

	lgm,s2 = np.loadtxt("lmsPkinputspec_new_ics.txt",unpack=True,usecols=[0,1])
	omegaz,deltac = np.loadtxt("flat.txt",unpack=True,usecols=[0,2])

	s2nu = np.interp(mass,lgm,s2)
	
	omega = cosmology.Planck13.Om(z)
	dc = np.interp(omega,omegaz,deltac)

	omegal = 1-omega
	gz = 5./2.*omega/(omega**(4./7.)-omegal+(1+omega/2.)*(1+omegal/70.))
	Dt = gz/(1+z)
	dc= dc/Dt

	nu = dc**2/s2nu
	nufnu = AA*(1.+1./(a*nu)**p)*np.sqrt(a*nu/(2*np.pi))*np.exp(-a*nu/2.)
	lgnu = np.log(nu)
        
	par = np.polyfit(mass,lgnu,2)
	#print par
	yyfit = par[2]+par[1]*mass+par[0]*mass**2
	der = np.zeros(nu.size)
	for j in range(1,nu.size):
	    der[j] = (yyfit[j]-yyfit[j-1])/(mass[j]-mass[j-1])

	mnz = 0.3*2.77e11*nufnu*der/2.3/10**mass
 	mnz_w= mnz*(1+m_c/10**mass)**(-1.3)

 	return mnz_w

def los_volume(z0,z2): # return the los total mass and tot_n halo in effective volume between (z0,z2)
	
	mass = np.linspace(6,9,200) # log10 mass range

	htot_f = cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h
	htot_b = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(zl))*cosmology.Planck13.h
	#print htot_f,htot_b

	zi = np.linspace(z0,z2,1000)
	dist = cosmology.Planck13.comoving_distance(zi)*cosmology.Planck13.h
	dist_mid = (dist[0]+dist[-1])/2
	#print dist_mid.value
	z_mid = np.interp(dist_mid.value,np.array(dist),zi) # redshift at the half diatance

	mnz = massf(z_mid,mass)

 	I_nlos = trapz(mnz,mass)
 	#I_mlos = trapz(mnz*10**mass,mass)

    ## then we calculate n_los on each z slice

	if np.logical_and(z0<zl,z2<zl): # foreground
		hstep1 = cosmology.Planck13.comoving_distance(z0)*cosmology.Planck13.h
		hstep2 = cosmology.Planck13.comoving_distance(z2)*cosmology.Planck13.h
		bstep1 = Re*hstep1/htot_f*f_hi
		bstep2 = Re*hstep2/htot_f*f_hi
		vol1 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(hstep2-hstep1)/3. ## f_hi vs f_low
		bstep1 = Re*hstep1/htot_f*f_low
		bstep2 = Re*hstep2/htot_f*f_low
		vol2 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(hstep2-hstep1)/3. 

		vol = vol1-vol2
		vol = vol.value

	elif np.logical_and(z0>zl,z2>zl): # background
		hstep1 = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(z0))*cosmology.Planck13.h
		hstep2 = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(z2))*cosmology.Planck13.h
		bstep1 = Re*hstep1/htot_b*f_hi
		bstep2 = Re*hstep2/htot_b*f_hi
		vol1 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(hstep1-hstep2)/3.
		bstep1 = Re*hstep1/htot_b*f_low
		bstep2 = Re*hstep2/htot_b*f_low
		vol2 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(hstep1-hstep2)/3. 

		vol = vol1-vol2
		vol=vol.value

	else: # zl is in the middle of z0,z2
		hstep1 = cosmology.Planck13.comoving_distance(z0)*cosmology.Planck13.h	
		bstep1 = Re*hstep1/htot_f*f_hi
		bstep2 = Re*f_hi
		vol1 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(htot_f-hstep1)/3.
		bstep1 = Re*hstep1/htot_f*f_low
		bstep2 = Re*f_low
		vol2 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(htot_f-hstep1)/3. 
		vol = vol1-vol2
	
		hstep2 = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(z2))*cosmology.Planck13.h
		bstep1 = Re*f_hi
		bstep2 = Re*hstep2/htot_b*f_hi
		vol1 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(htot_b-hstep2)/3.
		bstep1 = Re*f_low
		bstep2 = Re*hstep2/htot_b*f_low
		vol2 = np.pi*(bstep1**2 +bstep1*bstep2+bstep2**2)*(htot_b-hstep2)/3. 
		vol = vol + (vol1-vol2)
		vol=vol.value
		

	n_los = I_nlos*vol
	#m_los = I_mlos*vol # log10
	#m_los = 10**m_los

	return n_los

z = np.logspace(np.log10(0.2),np.log10(zs-0.2),200)
#print z
#m_zslide = np.zeros(len(z)-1)
n_zslide = np.zeros(len(z)-1)
zslide = np.zeros(len(z)-1)

for i in range(len(z)-1):
	n = los_volume(z[i],z[i+1])
	#m_zslide[i] = m
	n_zslide[i] = n
	zslide[i] = (z[i]+z[i+1])/2
	print (z[i]+z[i+1])/2,n

n_halo = int(np.sum(n_zslide))  # total LOS object

#print n_halo

#### start to draw LOS object (Here's the realization loop)

htot_f = cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h
htot_b = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(zl))*cosmology.Planck13.h
re_f = np.zeros(len(zslide))

mass = np.linspace(6,9,200)
mnz = np.zeros((len(zslide),200))
print "calculating mass function..."
for i in range(len(zslide)):
	
	mnz[i,:] = massf(zslide[i],mass) ## here

	# position
	if (zslide[i]<zl): # foreground
		hstep = cosmology.Planck13.comoving_distance(zslide[i])*cosmology.Planck13.h
		re_f[i] = Re_arc*hstep/htot_f

	else: # background
		hstep = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(zslide[i]))*cosmology.Planck13.h
		re_f[i] = Re_arc*hstep/htot_b

# draw z first
for j in range(n_real):

	print '####'
	print 'realization '+str(j+n_start)
	print '####'

	z_halo = np.random.choice(zslide,n_halo,p=n_zslide/np.sum(n_zslide))
	z_halo = np.sort(z_halo)

	n_draw = np.zeros(len(zslide))
	for i in range(len(zslide)):
		mask = (z_halo == zslide[i])
		n_draw[i] = len(z_halo[mask])

	n_draw = n_draw.astype(int)
	#print n_draw
	# then draw mass & x,y position (arcsec)
	
	m_halo = np.empty(0)
	x_list,y_list = np.empty(0),np.empty(0)
	for i in range(len(zslide)):
		if n_draw[i]>0:
			
			m_draw = np.random.choice(mass,n_draw[i],p=mnz[i,:]/np.sum(mnz[i,:]))
			m_halo = np.append(m_halo,m_draw)

			

			k=0
			while k<(n_draw[i]):
				#print "####"
				#print i
				#print '####'
				x_temp = np.random.rand(1)*2.0*f_hi*re_f[i]-f_hi*re_f[i]
				y_temp = np.random.rand(1)*2.0*f_hi*re_f[i]-f_hi*re_f[i]
				dis2_temp = np.sqrt(x_temp**2+y_temp**2)

				if np.logical_and(dis2_temp>(f_low*re_f[i]),(dis2_temp<f_hi*re_f[i])) :
					#print dis2_temp, x_temp,y_temp
					x_list = np.append(x_list,x_temp)
					y_list = np.append(y_list,y_temp)
					k=k+1
			##   position and mass drawing done

	m_halo = m_halo.flatten()
	m_halo = 10**m_halo

	r200_list = NFW.r200(z_halo,m_halo)
	c200_list = NFW.c200(z_halo,m_halo)

	print len(z_halo),len(m_halo),len(x_list),len(r200_list)


	real_i = path+real_folder+'los'+str(j+n_start)+'.txt'
	np.savetxt(real_i,np.c_[z_halo,m_halo,x_list,y_list,r200_list,c200_list])



