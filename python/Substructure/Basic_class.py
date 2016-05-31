""" Lens class & Cosmology class """
import numpy as np
import DistanceTool as distance

class Cosmology:
	def __init__(self,OmegaM=0.27,h0=0.704):
		self.OM = OmegaM
		self.h = h0

class Lens:
	def __init__(self,setup = np.nan):

		if setup == np.nan:
			print "No input lens setup file. Use default values."
			self.zl = 0.6
			self.zs = 2.0
			self.q = 0.2
			self.b = 0.5
			self.xc = 0.0
			self.yc = 0.0
			self.src_x = self.xc
			self.src_y = self.yc

		else:
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

	def obsData(self,data = np.NaN):
		if data == np.NaN:
			print 'Please select a obs_data file'

		else:
			self.img_x = data[:,0]
			self.img_y = data[:,1]
			self.img_f = data[:,2]
			self.img_err = data[:,3]
			self.img_ferr = data[:,4]

	def  lens_setFile(self,setup = np.NaN):

			self.zl = setup[0]
			self.zs = setup[1]
			self.q = setup[2]
			self.b = setup[3]
			self.xc = setup[4]
			self.yc = setup[5]
			self.src_x = setup[6]
			self.src_y = setup[7]
