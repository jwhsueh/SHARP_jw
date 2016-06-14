import numpy as np

### This code combine Lens properties of three projection

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

### Dandan's morphology pick
DanID_x = np.loadtxt('JenWei_table_disk_in_x.dat',dtype = 'int',unpack=True, usecols=[0])
DanID_y = np.loadtxt('JenWei_table_disk_in_y.dat',dtype = 'int',unpack=True, usecols=[0])
DanID_z = np.loadtxt('JenWei_table_disk_in_z.dat',dtype = 'int',unpack=True, usecols=[0])
DanID = np.union1d(DanID_x,DanID_y)
DanID = np.union1d(DanID,DanID_z)
#DanID.sort()

# If they are in DanID, morphology pick flag = 1

## Kinematics

catalog = basePath+'/kinematics_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])

bulge_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
disk_frac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'.dat'
Galaxy_ms = np.loadtxt(catalog2,dtype = 'float',unpack=True, usecols=[4])

## Lens properties from Dandan's full catalog [all mass range]

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_x.dat'
LensID = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
#Mass = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[1])/0.704

Re_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_x = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_y.dat'
#LensID_y = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Re_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_y = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

catalog3 = basePath+'/Dandan_Lens'+str(ssNumber)+'_z.dat'
#LensID_z = np.loadtxt(catalog3,dtype = 'int',unpack=True, usecols=[0])
Re_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[3])
DMfrac_z = np.loadtxt(catalog3,dtype='float',unpack=True,usecols=[4])

LensID = np.union1d(LensID_x,LensID_y)
LensID = np.union1d(LensID,LensID_z)

## get Re & DMfrac table
Re = np.zeros((LensID.size,3))*np.nan
DMfrac = np.zeros((LensID.size,3))*np.nan

for i in range(LensID.size):
	if LensID[i] in LensID_x:
		idx = list(LensID_x).index(LensID[i])
		Re[i,0] = Re_x[idx]
		DMfrac[i,0] = DMfrac_x[idx]

	if LensID[i] in LensID_y:
		idx = list(LensID_y).index(LensID[i])
		Re[i,1] = Re_y[idx]
		DMfrac[i,1] = DMfrac_y[idx]

	if LensID[i] in LensID_z:
		idx = list(LensID_z).index(LensID[i])
		Re[i,2] = Re_z[idx]
		DMfrac[i,2] = DMfrac_z[idx]


## Inclination angle

catalog4 = basePath+'/GalaxyInclination_'+str(ssNumber)+'.dat'
theta_x = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[1])
theta_y = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[2])
theta_z = np.loadtxt(catalog4,dtype='float',unpack=True,usecols=[3])

# 60 - 120 degree: high inclination, flag = 1

## The standard catalog use LensID

s_cata = open(basePath+'/Galaxy_Lens'+str(ssNumber)+'.dat','w')
s_cata.write('#[0]: SubfindID \n')
s_cata.write('#[1]: Mass [M_sun] \n')
s_cata.write('#[2]: Disk str fraction > 0.4, flag = 1 \n')
s_cata.write('#[3]: Bulge str fraction < 0.6, flag = 1 \n')
s_cata.write('#[4-6]: high inclination angle (60-120 degree), flag = 1; x,y,z \n')
s_cata.write('#[7-9]: Re x,y,z \n')
s_cata.write('#[10-12]: DM frac w/i Re x,y,z \n')
s_cata.write('#[13]: Dandan morphology pick, flag =1 \n')


for i in range(LensID.size):
	if LensID[i] in GalaxyID:
		idx = list(GalaxyID).index(LensID[i])

		## write everything about Galaxy catalog
		s_cata.write(str(LensID[i])+'	')
		s_cata.write(str(Galaxy_ms[idx])+'	')
		
		if disk_frac[idx]>0.4: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if bulge_frac[idx]<0.6: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if theta_x[idx]>60 and theta_x[idx]<120: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if theta_y[idx]>60 and theta_y[idx]<120: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		if theta_z[idx]>60 and theta_z[idx]<120: 
			s_cata.write('1  ')
		else: 
			s_cata.write('0  ')

		## Lens properties

		s_cata.write(str(Re[i,:])[1:-1]+'	')
		s_cata.write(str(DMfrac[i,:])[1:-1]+'	')

		## Dandan's morphology pick flag

		if LensID[i] in DanID:
			s_cata.write('1 \n')
		else:
			s_cata.write('0 \n')