import numpy as np
import scipy.stats as st
from astropy import cosmology
from scipy.integrate import trapz
import NFWprofile as NFW

def massf(z,mass,m_c):
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

def los_volume(z0,z2,light_cone,mass_info): # return the los total mass and tot_n halo in effective volume between (z0,z2)
	
	#mass = np.linspace(6,9,200) # log10 mass range
	zl,zs,Re,f_low,f_hi = light_cone[0],light_cone[1],light_cone[2],light_cone[3],light_cone[4]
	mass = mass_info[1]
	m_c = mass_info[0]

	Re = np.radians(Re_arc/3600.)*cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h

	htot_f = cosmology.Planck13.comoving_distance(zl)*cosmology.Planck13.h
	htot_b = (cosmology.Planck13.comoving_distance(zs)-cosmology.Planck13.comoving_distance(zl))*cosmology.Planck13.h
	#print htot_f,htot_b

	zi = np.linspace(z0,z2,1000)
	dist = cosmology.Planck13.comoving_distance(zi)*cosmology.Planck13.h
	dist_mid = (dist[0]+dist[-1])/2
	#print dist_mid.value
	z_mid = np.interp(dist_mid.value,np.array(dist),zi) # redshift at the half diatance

	mnz = massf(z_mid,mass,m_c)

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


