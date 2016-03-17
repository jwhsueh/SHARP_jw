
import numpy as np
import DistanceTool as distance
import scipy.integrate 

"""This function compute scaling radius of a NFW halo from velocity dispersion"""
""" Output is in arcsec """

def scaleR(cosmopara,lenspara,M_200,c_200):

## parameters
	q,b = lenspara.q, lenspara.b

	zl,zs =lenspara.zl, lenspara.zs 

	c = 3e8 # m/s
	G = 6.67e-11 # m^3/kg/s^2

	## velocity dispersion

	#Ds = distance.angular_distance(cosmopara, zs)
	Dl = distance.angular_distance(cosmopara, zl)
	#Dls = Ds - Dl

	## scale length

	Hz = 100*cosmopara.h*distance.Ez(cosmopara,zl) # km/s/Mpc
	Hz = Hz*1e3/3.08e22 # s-1

	rho_c = 3*Hz**2/(8*np.pi*G) # kg/m^3

	rs = (3*M_200/800./np.pi/rho_c/c_200**3)**(1.0/3.0) # m 
	rs_m = rs # m

	#rs = rs/3.08e22/cosmopara.h/100 # convert to Mpc
	#rs = np.degrees(rs/Dl)*3600 # arcsec 
	rs = distance.meter2arcs(cosmopara,rs_m,zl)

	return np.array([rs,rs_m])

""" velocity dispersion of lens [m/s] """

def velocity_dispersion(cosmopara,lenspara):

	c = 3e8 # m/s
	G = 6.67e-11 # m^3/kg/s^2

	zs = lenspara.zs
	zl = lenspara.zl

	Ds = distance.angular_distance(cosmopara, zs)
	Dl = distance.angular_distance(cosmopara, zl)
	Dls = Ds - Dl

	b,q = lenspara.b, lenspara.q
	b = np.radians(b/3600) # rad

	sig = b * np.sqrt(q)/(4*np.pi)*(Ds/Dls)*c**2 # velocity dispersion
	sig = np.sqrt(sig)

	return sig


def M200(cosmopara,lenspara,sig):

	c = 3e8 # m/s
	G = 6.67e-11 # m^3/kg/s^2

	zl = lenspara.zl

	Hz = 100*cosmopara.h*distance.Ez(cosmopara,zl) # km/s/Mpc
	Hz = Hz*1e3/3.08e22 # s-1

	rho_c = 3*Hz**2/(8*np.pi*G) # kg/m^3

	M_200 = sig**3*np.sqrt(3./800./np.pi/rho_c)*(2./G)**(3./2.)

	return M_200

def c200(h,zl,M_200):

	#M_200_sun = M_200/2e30 # M_sun

	c_200 = 5.71*(M_200/2e12*h)**(-0.084)*(1+zl)**(-0.47) # no unit

	return c_200



""" NFW P(r)*4 pi r^2, PDF """

def pdf(r,rs):
	x = r/rs
	return (4.0*np.pi*r**2)*1.0/x/(1+x)**2

""" Discrete CDF, ready for interpolate """

def cdf_d(ri,rs):

	#ri = np.linspace(0,r_end,r_end*10000)

	pdf_d = pdf(ri,rs) # discrete pdf
	cdf_d = np.zeros(len(ri))

	i = 1
	while i < len(ri):
		cdf_d[i] = cdf_d[i-1]+pdf_d[i]*(ri[1]-ri[0])
		i = i+1

	# normalization
	cdf_d = cdf_d/max(cdf_d)

	return cdf_d

def cdf(ri,rs):

	profile = lambda r: (4.0*np.pi*r**2)*1.0/(r/rs)/(1.+r/rs)**2

	integrate = np.zeros(len(ri))

	for i in range(len(ri)):

		I = scipy.integrate.quad(profile,0,ri[i])
		integrate[i] = I[0]

	return integrate

""" Draw from inverse_cdf will get cloned distribution """

def inverse_cdf(r,rs,r_end):

	ri = np.linspace(0,r_end,100)
	#Ix = cdf_d(ri,rs)
	Ix = cdf(ri,rs)
	Ix = Ix/Ix[-1]
	Iy = ri

	return np.interp(r,Ix,Iy)


""" Set up a class call halopara """
def set_halopara(cosmopara,lenspara):

	zl,zs = lenspara.zl, lenspara.zs
	
	''' expand a parameter to array contains 2 forms of units '''
	class halopara:
		sig = velocity_dispersion(cosmopara,lenspara) # m/s
		
		M_200 = np.array([M200(cosmopara,lenspara,sig),0.]) 
		M_200 = np.array([M_200[0],M_200[0]/2e30])# [kg, M_sun]
		
		c_200 = c200(cosmopara.h,lenspara.zl,M_200[1]) # no unit

		rs = scaleR(cosmopara,lenspara,M_200[0],c_200) # [arc sec, m]

		r_200 = rs*c_200 # [arc sec,m]

		rho_s = M_200[0]/(4.0*np.pi*rs[1]**3*(np.log(1.+c_200)-c_200/(1.+c_200))) # kg/m^3

		
	return halopara


""" Enclose mass within a radius ri [M_sun]"""
""" Input ri is in arc sec """

def enclose_mass(cosmopara,ri,zl,halopara):

	ri = distance.arcs2meter(cosmopara,ri,zl)

	rs = halopara.rs[1]  # m
	rho_s = halopara.rho_s # kg/m^3

	profile = lambda r: (4.0*np.pi*r**2)*1.0/(r/rs)/(1.+r/rs)**2

	M_en = np.zeros(len(ri))

	for i in range(len(ri)):
	
		I = scipy.integrate.quad(profile,0,ri[i])
		M_en[i] = rho_s*I[0]
	

	M_en = M_en/2e30 # M_sun

	#print 'enclose_mass'
	#print M_en

	return M_en

