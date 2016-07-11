import numpy as np
import DistanceTool as distance

basePath = '../../data/illustris_1/'

ssNumber = 99
redshift = 0.6

cataName = basePath+'Dandan_Phot'+str(ssNumber)+'_x.dat'
sigma_x = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[5])
R_x = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[4])

cataName = basePath+'Dandan_Phot'+str(ssNumber)+'_y.dat'
sigma_y = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[5])
R_y = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[4])

cataName = basePath+'Dandan_Phot'+str(ssNumber)+'_z.dat'
sigma_z = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[5])
R_z = np.loadtxt(cataName,dtype = 'float',unpack=True, usecols=[4])



''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27


DL = distance.luminosity_distance(cosmopara,redshift)*1e6 # pc

sigma_x = sigma_x+5.*(np.log10(DL)-1)
sigma_y = sigma_y+5.*(np.log10(DL)-1)
sigma_z = sigma_z+5.*(np.log10(DL)-1)

Sigma_x = 10**(-1.*sigma_x/2.5)
fx = Sigma_x/(np.pi*R_x**2)
mx = -2.5*np.log10(fx)

Sigma_y= 10**(-1.*sigma_y/2.5)
fy = Sigma_y/(np.pi*R_y**2)
my = -2.5*np.log10(fy)

Sigma_z = 10**(-1.*sigma_z/2.5)
fz = Sigma_x/(np.pi*R_z**2)
mz = -2.5*np.log10(fz)

disk = open(basePath+'DiskMag_'+str(ssNumber)+'.dat','w')

#Sigma0_x = Sigma_x*np.exp(1.)/2./(np.exp(1.)-2.)

#sigma0_x = -2.5*np.log10(Sigma0_x)
#print sigma0_x

#mu_x = sigma0_x+1.086*2
#print mu_x[mu_x>25.5].size
