""" Lens class & Cosmology class """
import numpy as np
import DistanceTool as distance

class Cosmology:
	def __init__(self,OmegaM=0.27,h0=0.704):
		self.OM = OmegaM
		self.h = h0

class Lens:
	def __init__(self,setup):

		self.lens_setFile(setup)


	def critical_density(self):
		c = 3e8 # m/s
		G = 6.67e-11 # m^3/kg/s^2

		## m to Mpc
		c = c/3.08e22  # Mpc/s
		G = G/(3.08e22)**3  # Mpc^3/kg/s^2

		cospara = Cosmology()

		Ds = distance.angular_distance(cospara,self.zs)
		Dl = distance.angular_distance(cospara,self.zl)
		Dls = Ds - Dl

		Sigma_c = c**2/(4.0*np.pi*G)*Ds/(Dl*Dls) # kg/Mpc^2
		Sigma_c = Sigma_c/2e30*cospara.h**2 # h M_sun/Mpc^2

		return Sigma_c

	def obsData(self,data):

		self.img_x = data[:,0]
		self.img_y = data[:,1]
		self.img_f = data[:,2]
		self.img_err = data[:,3]
		self.img_ferr = data[:,4]

	def  lens_setFile(self,setup):

			self.zl = setup[0]
			self.zs = setup[1]
			self.b = setup[2]
			self.xc = setup[3]
			self.yc = setup[4]
			self.q = setup[5]
			self.PA = setup[6]
			self.gamma1 = setup[7]
			self.gamma2 = setup[8]
			self.src_x = setup[9]
			self.src_y = setup[10]
