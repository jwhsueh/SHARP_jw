import numpy as np
import h5py
import DistanceTool as distance
import snapshot

### This code combine Lens properties of three projection

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

## get header file
headerPath = '/Volumes/narsil_1/jwhsueh/illustris_1'

f = h5py.File(snapshot.snapPath(headerPath,ssNumber),'r')
header = dict( f['Header'].attrs.items() )

redshift = header['Redshift']

### Dandan's morphology pick
DanID_x = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'int',unpack=True, usecols=[0])
DanID_y = np.loadtxt('JenWei_table_disk_in_y.dat',dtype = 'int',unpack=True, usecols=[0])
DanID_z = np.loadtxt('JenWei_table_disk_in_z.dat',dtype = 'int',unpack=True, usecols=[0])
DanID = np.union1d(DanID_x,DanID_y)
DanID = np.union1d(DanID,DanID_z)
#DanID.sort()

# If they are in DanID, morphology pick flag = 1

## Kinematics

catalog = basePath+'/kinematics_'+str(ssNumber)+'_sig.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

bulge_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
disk_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'_sig.dat'
GalaxyID = np.loadtxt(catalog2,dtype = 'int',unpack=True, usecols=[0])
Galaxy_ms = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[4])
Galaxy_str = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[5])
Galaxy_sig = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[6])
Galaxy_relax = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[8])
Galaxy_V = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[9])

## Lens properties from Dandan's full catalog [all mass range]

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_x.dat'
LensID_x = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
#Mass = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[1])/0.704

Re_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_y.dat'
LensID_y = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Re_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_z.dat'
LensID_z = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Re_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

LensID = np.union1d(LensID_x,LensID_y)
LensID = np.union1d(LensID,LensID_z)

## Photometry properties (Dandan's full catalog)
catalog3 = basePath+'/Dandan_Phot'+str(ssNumber)+'_x.dat'
Reff_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])
bright_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[5])

catalog3 = basePath+'/Dandan_Phot'+str(ssNumber)+'_y.dat'
Reff_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])
bright_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[5])

catalog3 = basePath+'/Dandan_Phot'+str(ssNumber)+'_z.dat'
Reff_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])
bright_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[5])


## absolute mag to apparent mag

''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

DL = distance.luminosity_distance(cosmopara,redshift)*1e6 # pc

bright_x = bright_x+5.*(np.log10(DL)-1) 
bright_y = bright_y+5.*(np.log10(DL)-1) 
bright_z = bright_z+5.*(np.log10(DL)-1) 


Sigma_x = 10**(-1.*bright_x/2.5)
fx = Sigma_x/(np.pi*Reff_x**2)
mx = -2.5*np.log10(fx)

Sigma_y= 10**(-1.*bright_y/2.5)
fy = Sigma_y/(np.pi*Reff_y**2)
my = -2.5*np.log10(fy)

Sigma_z = 10**(-1.*bright_z/2.5)
fz = Sigma_z/(np.pi*Reff_z**2)
mz = -2.5*np.log10(fz)

## get Re & DMfrac table
Re = np.zeros((LensID.size,3))*np.nan
DMfrac = np.zeros((LensID.size,3))*np.nan
bright = np.zeros((LensID.size,3))*np.nan

for i in range(LensID.size):
	if LensID[i] in LensID_x:
		idx = list(LensID_x).index(LensID[i])
		Re[i,0] = Re_x[idx]
		DMfrac[i,0] = DMfrac_x[idx]
		bright[i,0] = mx[idx]

	if LensID[i] in LensID_y:
		idx = list(LensID_y).index(LensID[i])
		Re[i,1] = Re_y[idx]
		DMfrac[i,1] = DMfrac_y[idx]
		bright[i,1] = my[idx]

	if LensID[i] in LensID_z:
		idx = list(LensID_z).index(LensID[i])
		Re[i,2] = Re_z[idx]
		DMfrac[i,2] = DMfrac_z[idx]
		bright[i,2] = mz[idx]


## Inclination angle

catalog4 = basePath+'/Inclination_'+str(ssNumber)+'_sig.dat'
#thetaID = np.loadtxt(catalog4,dtype = 'int',unpack=True, usecols=[0])
theta_x = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[1])
theta_y = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[2])
theta_z = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[3])

# 60 - 120 degree: high inclination, flag = 1

## The standard catalog use LensID

s_cata = open(basePath+'/Galaxy_Lens'+str(ssNumber)+'_sig.dat','w')
s_cata.write('#[0]: SubfindID \n')
s_cata.write('#[1]: Total Mass [M_sun] \n')
s_cata.write('#[2]: Stellar Mass [M_sun] \n')
s_cata.write('#[3]: Velocity Dispertion [km/s] \n')
# add CM here
s_cata.write('#[4]: Disk str fraction > 0.4, flag = 1 \n')
s_cata.write('#[5]: Bulge str fraction < 0.6, flag = 1 \n')
s_cata.write('#[6-8]: high inclination angle (60-120 degree), flag = 1; x,y,z \n')
s_cata.write('#[9-11]: Re x,y,z \n')
s_cata.write('#[12-14]: DM frac w/i Re x,y,z \n')
s_cata.write('#[15-17]: DIsk magnitude (r-band)  \n')
s_cata.write('#[18]: Dandan morphology pick, flag =1 \n')
s_cata.write('#[19]: relax = 0, not relax =1 \n')
s_cata.write('#[20]: Galaxy magnitude r-band\n')


for i in range(LensID.size):
	if LensID[i] in GalaxyID:
		idx = list(GalaxyID).index(LensID[i])

		## write everything about Galaxy catalog
		s_cata.write(str(LensID[i])+'	')
		s_cata.write(str(Galaxy_ms[idx])+'	')
		s_cata.write(str(Galaxy_str[idx])+'	')
		s_cata.write(str(Galaxy_sig[idx])+'	')
		
		if disk_frac[idx]>0.4: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if bulge_frac[idx]<0.6: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')
		

		if theta_x[idx]>80. and theta_x[idx]<100.: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if theta_y[idx]>80. and theta_y[idx]<100.: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if theta_z[idx]>80 and theta_z[idx]<100.: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		## Lens properties

		s_cata.write(str(Re[i,:])[1:-1]+'	')
		s_cata.write(str(DMfrac[i,:])[1:-1]+'	')
		s_cata.write(str(bright[i,:])[1:-1]+'	')

		## Dandan's morphology pick flag

		if LensID[i] in DanID:
			s_cata.write('1  ')
		else:
			s_cata.write('0  ')

		# relax flag
		s_cata.write(str(int(Galaxy_relax[idx]))+' ')

		# Photometry
		s_cata.write(str(Galaxy_V[idx])+'\n')