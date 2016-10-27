import numpy as np
import snapshot
import h5py
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt

#### --------------change path setting here -------------#####
## snapshotPath: where the snapshot directory is
## catalog: Galaxy catalog (read in)
## out_cata: Inclination catalog (write out)

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'
ssNum = 99

catalog='../../data/illustris_1/snapshot099_catalog/disk_099.dat'
out_cata = open('../../data/illustris_1/Rot_'+str(snapNum)+'_all.dat','w')

#### -------------------------------- #####

class cosmopara:
	h = 0.704
	OM = 0.27

## this function deal w/ galaxy on the boundary
def boundary(ci):
	
	#ci[ci< 0.5*boxsize] = ci[ci< 0.5*boxsize]+boxsize

	for i in ci:
		if i<0.5*boxsize:
			i = i+boxsize

	return ci

hubble_time=distance.Hubble_time(cosmopara) # Gyr
#### ---------------------------- ####

## read in snapshot header

f = h5py.File(snapshot.snapPath(snapshotPath,ssNum),'r')
header = dict( f['Header'].attrs.items() )

boxsize = header['BoxSize']  # ckpc/h
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

# age of universe
t_box=hubble_time*a**(3./2.)
print t_box

## read in Galaxy properties

subID=np.loadtxt(catalog)
subID=subID.astype(int)

pos=groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloPos']) # ckpc/h
pos=pos*a/cosmopara.h # kpc

# get peculiar velocity of group
SubVel = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloVel']) # km/s

R_half=groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloHalfmassRadType'])[:,4] # stellar half mass radius [ckpc/h]
R_half=R_half*a/cosmopara.h

# young star age cut
disk_age=5 # Gyr

Axis = np.zeros((len(subID),3))*np.nan  # principle rotational axis of each galaxy
theta_i = np.zeros((len(subID),3))*np.nan # inclination angle of three main axis

#subID=subID[2:3]
subID=np.array([6])

for i in range(len(subID)):

	subhalo_st=snapshot.loadSubhalo(snapshotPath,snapNum,subID[i],'stars')
	subhalo_gas=snapshot.loadSubhalo(snapshotPath,snapNum,subID[i],'gas')

	center=pos[subID[i],:]

	coord = subhalo_st['Coordinates']*a/cosmopara.h # kpc
	st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]
	coord = subhalo_gas['Coordinates']*a/cosmopara.h # kpc
	gas_x,gas_y,gas_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(st_x)-min(st_x)> 0.5*boxsize): 
		st_x = boundary(st_x)
		gas_x = boundary(gas_x)
	if (max(st_y)-min(st_y)> 0.5*boxsize): 
		st_y = boundary(st_y)
		gas_y = boundary(gas_y)
	if (max(st_z)-min(st_z)> 0.5*boxsize): 
		st_z = boundary(st_z)
		gas_z = boundary(gas_z)

	## change ref point to subhalo center
	st_pos=np.transpose(np.array([st_x-center[0],st_y-center[1],st_z-center[2]]))
	gas_pos=np.transpose(np.array([gas_x-center[0],gas_y-center[1],gas_z-center[2]]))

	## Hubble drag
	st_hd=st_pos*100.*distance.Ez(cosmopara,redshift)
	st_v = subhalo_st['Velocities']*np.sqrt(a) # km/s
	st_v=st_v+st_hd	# apply Hubble drag

	gas_hd=gas_pos*100.*distance.Ez(cosmopara,redshift)
	gas_v = subhalo_gas['Velocities']*np.sqrt(a) # km/s
	gas_v=gas_v+gas_hd

	st_ms = subhalo_st['Masses']*1e10/cosmopara.h
	gas_ms = subhalo_gas['Masses']*1e10/cosmopara.h

	# age of star
	#age_scale=subhalo_st['GFM_StellarFormationTime'] # in scale factor
	#st_age=t_box-hubble_time*age_scale**(3./2.) # Gyr

	st_U = subhalo_st['Potential']/a # (km/s)^2
	gas_U = subhalo_gas['Potential']/a # (km/s)^2

	## only use particles w/i 2x half mass radius
	st_r = np.sqrt(st_pos[:,0]**2+st_pos[:,1]**2+st_pos[:,2]**2)
	gas_r = np.sqrt(gas_pos[:,0]**2+gas_pos[:,1]**2+gas_pos[:,2]**2)
	r_cut = R_half[subID[i]]

	# apply raidus mask
	mask = st_r < 2.0*r_cut
	st_pos = st_pos[mask,:]
	st_v = st_v[mask,:]
	st_ms = st_ms[mask]
	st_U=st_U[mask]
	st_r=st_r[mask]

	mask = gas_r < 2.0*r_cut
	gas_pos = gas_pos[mask,:]
	gas_v = gas_v[mask,:]
	gas_ms = gas_ms[mask]
	gas_U=gas_U[mask]
	gas_r=gas_r[mask]

	'''
	## only use young star
	mask=st_age<disk_age
	st_x,st_y,st_z = st_x[mask],st_y[mask],st_z[mask]
	st_vx,st_vy,st_vz = st_vx[mask],st_vy[mask],st_vz[mask]
	st_ms = st_ms[mask]
	U=U[mask]
	st_age=st_age[mask]
	#print st_age[:100]
	'''

	st_L=np.cross(st_pos,st_v)
	gas_L=np.cross(gas_pos,gas_v)

	# use star particle to calculate principal axis
	L_len = np.sqrt(np.sum(st_L[:,0]*st_ms)**2+np.sum(st_L[:,1]*st_ms)**2+np.sum(st_L[:,2]*st_ms)**2)
	
	# angular momentum unit vector
	L_u=np.array([np.sum(st_L[:,0])/L_len,np.sum(st_L[:,1])/L_len,np.sum(st_L[:,2])/L_len])

	#Axis[i,:] = np.array([L_xu,L_yu,L_zu])

	theta=np.arccos(L_u)
	theta_i[i,:]=np.degrees(theta)

	print 'subID = '+str(subID[i])
	print 'Inclination angle = '+ str(theta_i[i,:])[1:-1]

	## project to pricipal axis [specific]
	st_Jz=np.dot(st_L,L_u)
	gas_Jz=np.dot(gas_L,L_u)

	#st_Jz=-1.0*(st_Lx*np.sin(theta_x+np.pi/2)+st_Ly*np.sin(theta_y+np.pi/2)+st_Lz*np.sin(theta_z+np.pi/2))/st_ms
	#gas_Jz=-1.0*(gas_Lx*np.sin(theta_x+np.pi/2)+gas_Ly*np.sin(theta_y+np.pi/2)+gas_Lz*np.sin(theta_z+np.pi/2))/gas_ms

	## specific binding energy
	K=0.5*(st_v[:,0]**2+st_v[:,1]**2+st_v[:,2]**2)
	# centrifugal potential
	U_cen=0.5*(st_L[:,0]**2+st_L[:,1]**2+st_L[:,2]**2)/(st_r**2)

	Es=st_U+K-U_cen
	#print Es
	#Es=gas_U

	mask=Es<0
	st_Jz=st_Jz[mask]
	Es=Es[mask]

	E_bin=np.linspace(np.min(Es)/1e5,np.max(Es)/1e5,1001)
	print np.min(Es),np.max(Es)
	E_list=np.zeros(len(E_bin)-1)
	Jz_list=np.zeros(len(E_bin)-1)


	for j in range(len(E_bin)-1):
		mask=Es/1e5<E_bin[j+1]

		Jz_pool=st_Jz[mask]/1e3
		#Jz_pool=gas_Jz[mask]/1e3
		Es_pool=Es[mask]/1e5

		Jz_binmax=np.max(Jz_pool)
		idx=list(Jz_pool).index(Jz_binmax)
		Es_binmax=Es_pool[idx]

		Jz_list[j]=Jz_binmax
		E_list[j]=Es_binmax

	## fit with quadratic curve

	coef=np.poly1d(np.polyfit(E_list,Jz_list,2))
	#cort=coef(np.min(Es)/1e5)
	cort=0
	
	plt.scatter(Es/1e5,st_Jz/1e3,marker='.',s=1)
	#plt.scatter(Es/1e5,gas_Jz/1e3,marker='.',s=1)
	#plt.scatter(E_list,Jz_list,c='r')
	plt.plot(E_bin,coef(E_bin)-cort,'-r')
	plt.xlabel("Binding Energy E/(10^5 km^2/s^2)")
	plt.ylabel("Jz/(10^3 kpc km/s)")
	plt.title('Star Particles')
	plt.show()
	

	## calculate circularity

	#J_circ=coef(st_U/1e5)-cort
	J_circ=coef(st_U/1e5)
	print J_circ
	eps=st_Jz/1e3/J_circ

	# histogram
	hist,bin_edge=np.histogram(eps,bins=100)
	idx=list(hist).index(np.max(hist))
	#eps_2=eps-bin_edge[idx]
	eps_2=eps

	plt.hist(eps_2,bins=30)
	#plt.show()

	## kinematics component star fraction
	thin_disk=eps_2[eps_2>0.7]
	bulge=eps_2[eps_2<0]
	#bulge=eps_2[eps_2<bin_edge[idx+1]]

	#print thin_disk
	#print bulge

	print float(len(thin_disk))/len(eps),2*float(len(bulge))/len(eps)
	