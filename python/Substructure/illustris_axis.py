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

subID=np.array([6])

for i in range(len(subID)):

	subhalo_st=snapshot.loadSubhalo(snapshotPath,snapNum,subID[i],'stars')

	center=pos[subID[i],:]

	coord = subhalo_st['Coordinates']*a/cosmopara.h # kpc
	st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

	if (max(st_x)-min(st_x)> 0.5*boxsize): 
		st_x = boundary(st_x)
	if (max(st_y)-min(st_y)> 0.5*boxsize): 
		st_y = boundary(st_y)
	if (max(st_z)-min(st_z)> 0.5*boxsize): 
		st_z = boundary(st_z)

	## change ref point to subhalo center
	st_x,st_y,st_z = st_x-center[0],st_y-center[1],st_z-center[2]

	## Hubble drag
	hd_x,hd_y,hd_z = st_x*100.*distance.Ez(cosmopara,redshift),st_y*100.*distance.Ez(cosmopara,redshift),st_z*100.*distance.Ez(cosmopara,redshift)

	vel = subhalo_st['Velocities']*np.sqrt(a) # km/s

	st_vx,st_vy,st_vz = vel[:,0],vel[:,1],vel[:,2]
	# apply Hubble drag
	st_vx,st_vy,st_vz = st_vx+hd_x,st_vy+hd_y,st_vz+hd_z

	st_ms = subhalo_st['Masses']*1e10/cosmopara.h

	# age of star
	age_scale=subhalo_st['GFM_StellarFormationTime'] # in scale factor
	st_age=t_box-hubble_time*age_scale**(3./2.) # Gyr

	U = subhalo_st['Potential']/a # (km/s)^2

	## only use particles w/i 2x half mass radius
	st_r = np.sqrt(st_x**2+st_y**2+st_z**2)
	r_cut = R_half[subID[i]]
	mask = st_r < 2.0*r_cut

	# apply raidus mask
	st_x,st_y,st_z = st_x[mask],st_y[mask],st_z[mask]
	st_vx,st_vy,st_vz = st_vx[mask],st_vy[mask],st_vz[mask]
	st_ms = st_ms[mask]
	st_age=st_age[mask]
	U=U[mask]

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

	st_Lx = (st_y*st_vz-st_z*st_vy)*st_ms
	st_Ly = (st_z*st_vx-st_x*st_vz)*st_ms
	st_Lz = (st_x*st_vy-st_y*st_vx)*st_ms

	L_len = np.sqrt(np.sum(st_Lx)**2+np.sum(st_Ly)**2+np.sum(st_Lz)**2)

	L_xu,L_yu,L_zu = np.sum(st_Lx)/L_len,np.sum(st_Ly)/L_len,np.sum(st_Lz)/L_len  # unit vector component

	Axis[i,:] = np.array([L_xu,L_yu,L_zu])

	theta_x,theta_y,theta_z = np.arccos(L_xu),np.arccos(L_yu),np.arccos(L_zu)

	theta_i[i,:] = np.array([np.degrees(theta_x),np.degrees(theta_y),np.degrees(theta_z)])

	print 'subID = '+str(subID[i])
	print 'Inclination angle = '+ str(theta_i[i,:])[1:-1]

	## project to pricipal axis [specific]
	st_Jz=(st_Lx*L_xu+st_Ly*L_yu+st_Lz*L_zu)/st_ms

	## specific binding energy
	st_uz=st_vx*np.sin(theta_x)+st_vy*np.sin(theta_y)+st_vz*np.sin(theta_z)
	#Es = U+0.5*(st_uz/10)**2
	Es=U

	E_bin=np.linspace(np.min(Es)/1e5,np.max(Es)/1e5,1001)
	print np.min(Es),np.max(Es)
	E_list=np.zeros(len(E_bin)-1)
	Jz_list=np.zeros(len(E_bin)-1)


	for j in range(len(E_bin)-1):
		mask=Es/1e5<E_bin[j+1]

		Jz_pool=st_Jz[mask]/1e3
		Es_pool=Es[mask]/1e5

		Jz_binmax=np.max(Jz_pool)
		idx=list(Jz_pool).index(Jz_binmax)
		Es_binmax=Es_pool[idx]

		Jz_list[j]=Jz_binmax
		E_list[j]=Es_binmax

	## fit with quadratic curve

	coef=np.poly1d(np.polyfit(E_list,Jz_list,2))
	
	#plt.scatter(U/1e5,st_Jz/1e9,marker='.',s=1)
	plt.scatter(Es/1e5,st_Jz/1e3,marker='.',s=1)
	#plt.scatter(E_list,Jz_list,c='r')
	plt.plot(E_bin,coef(E_bin),'-r')
	plt.xlabel("Binding Energy E/(10^5 km^2/s^2)")
	plt.ylabel("Jz/(10^3 kpc km/s)")
	plt.title('Star Particles')
	plt.show()
	

	## calculate circularity

	J_circ=coef(Es/1e5)
	print J_circ
	eps=st_Jz/1e3/J_circ
	#print eps[:100]
	plt.hist(eps)
	#plt.show()

	## kinematics component star fraction
	thin_disk=eps_2[eps_2>0.7]
	bulge=eps_2[eps_2<0]
	#bulge=eps_2[eps_2<bin_edge[idx+1]]

	#print thin_disk
	#print bulge

	print float(len(thin_disk))/len(eps),2*float(len(bulge))/len(eps)
